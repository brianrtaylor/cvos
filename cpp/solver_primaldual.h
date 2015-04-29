//! @file solver_primaldual.h 
//! Code in this file solves the optimization problem of the form:
//! \f$ \min_{c:c\geq 0} \tau^T c + \tau_2\|c\|_{\infty} + \|W Dc\|_1 + \lambda^T \max(0, 1-D_{occ} c) + \kappa^T \max(0, 1-c) \f$.
//! Solution is obtained using the primal-dual method of Pock and Chambolle.

#ifndef SOLVER_PRIMALDUAL_H_
#define SOLVER_PRIMALDUAL_H_

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>
#include <limits>
#include <cassert>
#include <cstring>

#include "utils.h"
#include "sparse_op.h"
#include "bcv_diff_ops.h"


#ifdef HAVE_SSE
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif

// ensures output is printed in MATLAB.
#ifdef HAVE_MATLAB
#include "mex.h"
#define printf mexPrintf
#endif


using namespace std;

//! primal-dual *optimization* parameters (substructure of p.d. parameters)
struct pd_opt_params {
    //! parameter in prox. operators. should be in [0, sqrt(8)/2] (see Pock and Chambolle)
    float sigma_y;
    //! parameter in prox. operators. should be in [0, sqrt(8)/2] (see Pock and Chambolle)
    float sigma_c;
    //! maximum number of iterations
    int max_iterations;
    //! threshold for stopping optimization; stops when the value of 
    //! \$ |f(x_{t+1})-f(x_t)|/f(x_{t+1}) \$ falls below this number
    float fx_tolerance;
    //! threshold for stopping optimization without evaluating f(x) (which is
    //! somewhat costly). solver stops when \$ |x_{t+1} - x_t |/x_{t+1} \$
    //! falls below this number
    float dx_tolerance;
    //! specified how often to output optimization progress.
    //! set to 0 for silent optimization.
    int verbosity;

    //! When set to 0, L1 regularization used for layer values (as in
    //! "Semantic segmentation exploiting occlusion relations within a convex
    //! optimization framework" EMMCVPR13).
    //! Otherwise, when set to 1, L-infinity regularization is used (as in
    //! "Detachable object detection", TPAMI2012).
    bool use_ell_infinity;

    // when set to 0, we look at 'use_ell_infinity' to decide whether to use L1
    // or L-infinity. otherwise, both penalties are used.
    bool use_both_penalties;

    //! Upper bound on the number of layers. When set to <= 0, it is not used.
    float layer_upper_bound;

    // when set to 1, the code expects the field 'kappa' to be of size nnodes x 1
    // that field specifies the 'foreground prior'. more precisely, it is an
    // additional term in the objective, of the form:
    //          kappa^T max(0, 1 - c )
    // (i.e. a term that forces c >= 1, or keeps c = 0 at a cost kappa.)
    // Note that from the perspective of optimization, this is enforced slightly
    // differently depending on whether L1 or L-INF penalty is used.
    // This objective is somwhat easier to enforce with L1
    bool use_temporal_penalty;

    //! initial point used in the optimization. if 'use_initial_layers' set to 1,
    //! then, 'init_layers' is used to hot-start the problem. otherwise, 
    //! if 'use_initial_layers'=0, then init_layers is not used.
    vector<float> init_layers;
    bool use_initial_layers;
#ifdef HAVE_MATLAB
    // if set to nonzero value, optimization will not be stopped until
    // abs(f(x)-fx_truth)/fx_truth < fx_truth_eps, or a very large number
    // of iterations has passed.
    float fx_truth;
    float fx_truth_eps;
#endif
}; 

//! primal-dual *problem* parameters (substructure of p.d. parameters)
struct pd_prob_params {
    //! weight on L1 regularization
    vector<float> tau1;
    //! weight on L-infinity regularization
    float tau2;
    
    //! vector weight on foreground. this is used ONLY IF the field
    //! 'use_temporal_penalty' is set to 1.
    vector<float> kappa;

    float layer_ub; //! layer upper bound
    
    int nnodes;
    int nedges;
    int nocc_constraints;

    int rows;
    int cols;
    sparse_op<int> D;
    sparse_op<int> Docc;
    vector<float> weights;
    vector<float> occweights;
    //! When set to 0, (D,Docc) are used (which represent generic sparse linear
    //! operators. When set to 1, (Dx,Dy, and Docc) are used. This can only be
    //! done if optimization proceeds on a pixel lattice.
    //! When possible, it is best to set this to 1.
    bool solve_pixelwise;
};

//! parameters for the primal-dual solver (main structure)
struct primaldual_params {
    pd_opt_params opt;
    pd_prob_params prob;
};

//! primal-dual solver of detachable object detection problem
class solver_primaldual {
public:
    primaldual_params* params;
    
    solver_primaldual();
    ~solver_primaldual();
    solver_primaldual(primaldual_params* p); 

    //! Set initial point
    void initialize(vector<float>& x);
   
    //! Returns zero value if parameters are bad, 1 if parameters are OK 
    int check_params(); 
    //! Solves optimization and returns a vector of layers.
    vector<float> solve(); 

    //! Evaluates the optimization problem function
    float eval_opt_func_value();
    //! Evaluates the optimization problem function (arguments are for precomputed data)
    float eval_opt_func_value( const vector<float>& Dx, const vector<float>& Doccx);

    void print_cost();

    //! estimates || D^T D ||_2 using the power method.
    float estimate_divgrad_norm();
private:
    float simplex_proj_lb;
    float simplex_proj_ub;
    float simplex_proj_sol;
    int simplex_num_iters;
    float simplex_scaling;
    vector<float> sol_c;
    vector<float> sol_y1;
    vector<float> sol_y2;
    vector<float> sol_y3;
    vector<float> sol_x;

    void show_warning();

    void preadjust_algweights();

    // projects the vector onto the probability simplex. 
    void project_onto_prob_simplex(vector<float>& x, float eta, float tol=1e-5f);
        
    // evaluates: sum( max(x-alpha, 0) ) - eta
    float eval_simplex_disagreement(const vector<float>& x, float alpha, float eta);
    
    float eval_simplex_disagreement(const vector<float>& x, const vector<float>& alpha, float eta);

    void solve_y1(vector<float>& DX, float sigma);
    void solve_y3(float sigma);
 
    // -----------------------------------------------------------------------------------------
    //                      functions with a vector for L1 regularization
    // -----------------------------------------------------------------------------------------
    float solve_c(vector<float>& DTY, vector<float>& DOCCTY, 
            float sigma, vector<float>& tau);    
    float solve_c_temporal(vector<float>& DTY, vector<float>& DOCCTY, 
            float sigma, vector<float>& tau);
    float solve_c_ell_infinity(vector<float>& DTY, vector<float>& DOCCTY, float sigma);
    float solve_c_ell1_ell_infinity(vector<float>& DTY, vector<float>& DOCCTY, float sigma);
    void ell1_ell_infinity_recover_neg(vector<float>& c_out, 
                                const vector<float>& c, float sigma, const vector<float>& tau);

    void to_nonnegative_orthant(vector<float>& x);
     
    #ifdef HAVE_SSE
    void solve_y1_sse(vector<float>& DX, float sigma);
    void solve_y3_sse(float sigma);
    // -----------------------------------------------------------------------------------------
    //                      functions with a vector for L1 regularization
    // -----------------------------------------------------------------------------------------
    float solve_c_sse(vector<float>& DTY, vector<float>& DOCCTY, 
            float sigma, vector<float>& tau);
    float solve_c_temporal_sse(vector<float>& DTY, vector<float>& DOCCTY, 
            float sigma, vector<float>& tau);
    // -----------------------------------------------------------------------------------------
    void solve_c_ell_infinity_prepare_(vector<float>& out, 
            vector<float>& DTY, vector<float>& DOCCTY, float sigma, bool temporal,
            int project_onto_feasible_set = 1); 
    float solve_c_ell_infinity_update_(
            vector<float>& c_temp, vector<float>& c_simplex);
    float eval_simplex_disagreement_sse(const vector<float>& x, float alpha, float eta);
    float eval_simplex_disagreement_sse(const vector<float>& x, const vector<float>& alpha, float eta);
    #endif
};

#endif
