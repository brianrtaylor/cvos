#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include "../cpp/solver_primaldual.h"
#include "../cpp/bcv_diff_ops.h"
#include "../cpp/sparse_op.h"
#include "../cpp/utils.h"

#define throw_error(x) mexPrintf("'%s' does not exist.\n", (x) ); return;
#define throw_error_double(x) mexPrintf("'%s' is not double.\n", (x) ); return;
#define MAT2C(x) ((x)-1)
#define C2MAT(x) ((x)+1)

// ------------------------------------------------------------------------
//                  USAGE: layers = pd_wrapper(problem)
// ------------------------------------------------------------------------
using namespace std;

void print_usage();
template<typename T> double mxGetScalarFromStruct(T* value, 
                        const mxArray* problem, const char* field_name);
template<typename T> void mxGetVectorFromStruct(vector<T>& x, 
                        const mxArray* problem, const char* field_name);
bool mxFieldExists(const mxArray* problem, const char* field_name);
primaldual_params init_params_struct();
void print_primaldual_params(const primaldual_params& params);
int is_double(const mxArray* pm) { return (mxDOUBLE_CLASS == mxGetClassID(pm)); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // --------------------------------------------------------------------
    if ( (nrhs != 1) || (nlhs > 2)) {
        mexPrintf("incorrect number of arguments.\n\n");
        print_usage();
        return;
    }
    const mxArray* problem = prhs[0];
    const mxArray* data;
    if (!mxIsStruct(problem)) { 
        mexPrintf("Input is not a problem structure.\n");
        return;
    }
    double *imsize;
    // default parameters
    int rows = 0;
    int cols = 0;
    // --------------------------------------------------------------------
    //              get image size (if data is on the lattice)
    // --------------------------------------------------------------------    
    data = mxGetField(problem, 0, "imsize");
    if (data != NULL) {
        if (!is_double(data)) { throw_error_double("imsize") };
        imsize = (double*)mxGetPr( data );
        rows = imsize[0];
        cols = imsize[1];
    }
    //---------------------------------------------------------------------
    //              Initialize default parameters
    //---------------------------------------------------------------------
    primaldual_params params = init_params_struct();
    
    // Make sure that several *required* fields are present.
    if (!mxFieldExists(problem, "Wx")) { return; }
    if (!mxFieldExists(problem, "lambda")) { return; }
    if (!mxFieldExists(problem, "constraints")) { return; }
    if (!mxFieldExists(problem, "edges")) { return; }
    if (!mxFieldExists(problem, "tau1")) { return; }
    
    mxGetVectorFromStruct(params.prob.weights, problem, "Wx");
    // ------------------------------------------------------------------//                    
    data = mxGetField(problem, 0, "Wx");
    if ( ((int)mxGetN(data) > 1) && ((int)mxGetM(data) > 1) ) {
        mexPrintf("'Wx' should be a vector, not a matrix.\n");
        return;
    }
    // ------------------------------------------------------------------//                    
    mxGetVectorFromStruct(params.prob.occweights, problem, "lambda");
    // ------------------------------------------------------------------//                    
    data = mxGetField(problem, 0, "lambda"); 
    if (((int)mxGetN(data) > 1) && ((int)mxGetM(data) > 1)) {
        mexPrintf("'lambda' should be a vector, not a matrix.\n");
        return;
    }
    // ------------------------------------------------------------------//                    
    mxGetVectorFromStruct(params.prob.tau1, problem, "tau1");    
    if (params.prob.tau1.size() != rows*cols) { 
        mexPrintf("error: tau1.size() = %d, rows=%d, cols=%d\n", 
                params.prob.tau1.size(), rows, cols);
        return;
    }    
    
    data = mxGetField(problem, 0, "constraints");
    int nocc_constraints = (int)mxGetM(data);    
    // ------------------------------------------------------------------//
    if (((int)mxGetM(data) == 2) || ((int)mxGetN(data) !=2)) { 
        mexPrintf("'constraints' should be a Nx2 array.\n");
        return;
    }
    vector<int> occluded = vector<int>(nocc_constraints);
    vector<int> occluder = vector<int>(nocc_constraints);
    double* constraints_ = (double*)mxGetPr(data);
    for (size_t i = 0; i < nocc_constraints; ++i) { 
        occluder[i] = MAT2C( (int)constraints_[i] );
        occluded[i] = MAT2C( (int)constraints_[i + nocc_constraints] );
    }
    // --------------------------------------------------------------------
    //                          get scalar fields
    // --------------------------------------------------------------------
    mxGetScalarFromStruct(&params.prob.nnodes, problem, "nnodes");
    mxGetScalarFromStruct(&params.prob.nedges, problem, "nedges");
    mxGetScalarFromStruct(&params.prob.tau2, problem, "tau2");
    mxGetScalarFromStruct(&params.prob.solve_pixelwise, problem, "SOLVE_PIXELWISE");

    mxGetScalarFromStruct(&params.opt.verbosity, problem, "verbosity");
    mxGetScalarFromStruct(&params.opt.fx_tolerance, problem, "fx_tolerance");
    mxGetScalarFromStruct(&params.opt.dx_tolerance, problem, "dx_tolerance");
    mxGetScalarFromStruct(&params.opt.max_iterations, problem, "max_iterations");
  
    mxGetScalarFromStruct(&params.opt.layer_upper_bound, problem, "layer_upper_bound");
    mxGetScalarFromStruct(&params.opt.use_both_penalties, problem, "use_both_penalties");
    mxGetScalarFromStruct(&params.opt.use_ell_infinity, problem, "USE_ELL_INFINITY");
 
    mxGetScalarFromStruct(&params.opt.sigma_c, problem, "sigma_c");
    mxGetScalarFromStruct(&params.opt.sigma_y, problem, "sigma_y");
    mxGetScalarFromStruct(&params.opt.fx_truth, problem, "fx_truth");
    mxGetScalarFromStruct(&params.opt.fx_truth_eps, problem, "fx_truth_eps");

    if (mxGetField(problem, 0, "init_layers")!=NULL) {
        int n = mxGetNumberOfElements( mxGetField(problem, 0, "imsize") );
        if (n == params.prob.nnodes) {
            params.opt.use_initial_layers = 1;
            mxGetVectorFromStruct(params.opt.init_layers, problem, "init_layers");
        } else {
            params.opt.use_initial_layers = 0;
        }
    }
    
    mxGetScalarFromStruct(&params.opt.use_temporal_penalty, problem, "USE_TEMPORAL_PENALTY");
    if (params.opt.use_temporal_penalty) { 
        if (!mxFieldExists(problem, "kappa")) {
            mexPrintf("USE_TEMPORAL_PENALTY=1, but kappa d.n.e.\n");
            throw_error("kappa");
        } else {
            mxGetVectorFromStruct(params.prob.kappa, problem, "kappa");
            if (params.prob.kappa.size() != params.prob.nnodes) {
               mexPrintf("kappa.size() does not match nnodes.\n");
               return;
            }
        }   
    }
    // --------------------------------------------------------------------
    //  should the problem should be solved on a graph or on the lattice. 
    // --------------------------------------------------------------------
    if ((params.prob.solve_pixelwise == 1) && ((rows == 0) || (cols==0))) {
        mexPrintf("SOLVE_PIXELWISE = 1 but rows=0 or cols=0.\n");
        mexPrintf("If you want to solve the problem on a *graph*, then\n");
        mexPrintf("specify SOLVE_PIXELWISE=0 and set 'nnodes' and 'nedges'\n");
        mexPrintf("If you want to solve the problem on a *lattice*, then\n");
        mexPrintf("set SOLVE_PIXELWISE=1, but specify 'rows' and 'cols'\n");
        mexPrintf("As it is right now, it is not clear how the problem should be solved.\n");
        return;
    }
    if ((params.prob.solve_pixelwise ==0) && 
                ((params.prob.nnodes ==0) || (params.prob.nedges==0)) ) {
        mexPrintf("SOLVE_PIXELWISE=0, but nnodes=0 or nedges=0.\n");
        mexPrintf("To solve the problem on a *graph*, you need to specify\n");
        mexPrintf("both 'nnodes' and 'nedges' in addition to 'SOLVE_PIXELWISE=1'\n");
        return;
    } 
    
    // --------------------------------------------------------------------
    if (params.prob.solve_pixelwise) { 
        mexPrintf("Solving problem on a pixel lattice.\n");
        params.prob.Docc = create_diff_op_from_data(occluded, occluder, rows*cols); 
        params.prob.nnodes = rows*cols;
        params.prob.nedges = rows*cols*2;
        params.prob.nocc_constraints = params.prob.Docc.nrows;
        params.prob.rows = rows;
        params.prob.cols = cols;
    } else {
        mexPrintf("Solving problem on a graph.\n");
        // this is only needed for superpixel graph type of problems.
        vector<int> edges1 = vector<int>(params.prob.nedges);
        vector<int> edges2 = vector<int>(params.prob.nedges);
        if (!mxFieldExists(problem, "edges")) { return; }
        double* edges = (double*)mxGetPr( mxGetField(problem, 0, "edges") );
        for (size_t i = 0; i < params.prob.nedges; ++i) { 
            edges1[i] = MAT2C( edges[i] );
            edges2[i] = MAT2C( edges[i+params.prob.nedges] );
        }        
        params.prob.Docc = create_diff_op_from_data(occluded, occluder, params.prob.nnodes); 
        params.prob.D = create_diff_op_from_data(edges1, edges2, params.prob.nnodes);
        params.prob.nnodes = params.prob.nnodes;
        params.prob.nedges = params.prob.nedges;
        params.prob.nocc_constraints = params.prob.Docc.nrows;
        params.prob.rows = 0;
        params.prob.cols = 0; 
    }    
    print_primaldual_params(params);    
    solver_primaldual pd(&params);
    
    vector<float> layers = vector<float>(params.prob.nnodes);    
    layers = pd.solve();

    pd.print_cost();
    // --------------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(params.prob.nnodes, 1, mxREAL);
    double* out = (double*)mxGetPr(plhs[0]);
    for (size_t i = 0; i < layers.size(); ++i) {
        out[i] = (double)layers[i];
    }
}

void print_usage() {
    mexPrintf("USAGE: layers = pd_wrapper(problem)\n");
    mexPrintf("  where problem is a structure containing (at least) the following fields:\n");
    mexPrintf("\n");
    mexPrintf("Wx - vector of weights defined on 'edges' (NEDGES x 1)\n");
    mexPrintf("edges - indices of edges (NEDGES x 2)\n");
    mexPrintf("lambda - vector of constraint weights (NOCC x 1)\n");
    mexPrintf("constraints - indices of constraints (NOCC x 2)\n");
    mexPrintf("nnodes - number of nodes in the problem\n");
    mexPrintf("nedges - number of edges in the problem\n");
    mexPrintf("imsize - (rows, cols)\n");
    mexPrintf("tau1 - weight on the L1 penalty (NNODES x 1)\n");
    mexPrintf("tau2 - weight on the Linf penalty (by default, unused)\n");
    mexPrintf("USE_ELL_INFINITY   - if '0', L1 norm is used. \n");
    mexPrintf("                   - if '1', Linf norm is used.\n");
    mexPrintf("USE_BOTH_PENALTIES - if '1', L1+Linf is used. \n");
    mexPrintf("                   - if '0', L1 or Linf is used (depending on other settings).\n");
    mexPrintf("layer_upper_bound - if set to 0, upper bound not used\n");
    mexPrintf("USE_TEMPORAL_PENALTY - 0/1 whether 'foreground prior is used.\n");
    mexPrintf("kappa - weight of the foreground prior. (NNODES x 1) \n");
    mexPrintf("\n");
    mexPrintf("max_iterations - maximum number of iterations\n");
    mexPrintf("fx_tolerance - function change tolerance (as stopping criterion)\n");
    mexPrintf("  when set to 0, function value is not calculated at all\n");
    mexPrintf("dx_tolerance - point change tolerance (as stopping criterion)\n");
    mexPrintf("verbosity - how often to print diagnostics\n");
    mexPrintf("init_layers - initial point in the optimization\n");
    mexPrintf("\n");
    mexPrintf("SOLVE_PIXELWISE - 0/1. If set to 0, problem is solved on a graph\n");
    mexPrintf("  if set to 1, it is assumed that the problem is specified on a pixel\n");
    mexPrintf("  lattice. In particular it is assumed that weights and edges are\n");
    mexPrintf("  in specified in a particular order: *edges = [dx_inds; dy_inds]*\n");
    mexPrintf("  this assumption is exploited to create fast difference operators.\n");
    
    mexPrintf("\n");
    mexPrintf("  Whenever possible (when the problem indeed corresponds to optimization\n");
    mexPrintf("  over a pixel lattice), you should be setting 'SOLVE_PIXELWISE=1'\n");
    mexPrintf("  In doing so, make sure that edges and Wx are specified correctly\n");
    mexPrintf("  (meaning that the assumption is satisfied).\n");    
    mexPrintf("\n");
}


//! Gets a scalar value from mx structure, and rewrites "T* value" with it.
//!
//! If 'field_name' exists in 'problem', 'value' is set to the first entry.
//! Otherwise, value is unchanged.
template<typename T> double mxGetScalarFromStruct(T* value, 
                        const mxArray* problem, const char* field_name) {
    const mxArray* data = mxGetField(problem, 0, field_name);
    if (data != NULL) {
        double* temp = (double*)mxGetPr( data );
        *value = (T)temp[0];
    }
}

//! Returns 0 if the field does not exist in the structure or is not double.
bool mxFieldExists(const mxArray* problem, const char* field_name) {
    const mxArray* data = mxGetField(problem, 0, field_name); 
    if (data == NULL) { 
        mexPrintf("'%s' does not exist.\n", field_name ); 
        return 0;
    } 
    if (!is_double(data)) { 
        mexPrintf("'%s' is not double.\n", field_name);
        return 0;
    }
    return 1;
}

//! Get vector from 'field_name' into 'x'.
template<typename T> void mxGetVectorFromStruct(vector<T>& x, 
                        const mxArray* problem, const char* field_name) {
    const mxArray* data = mxGetField(problem, 0, field_name); 
    if (data != NULL) {
        double* temp = (double*)mxGetPr(data);
        int n = mxGetNumberOfElements(data);
        vector<T> y( temp, temp + n);
        x = y;
    }
}


primaldual_params init_params_struct() {
    primaldual_params params;
    params.opt.sigma_y = 0.9f/sqrt(8.0f);
    params.opt.sigma_c = 0.9f/sqrt(8.0f);
    params.opt.max_iterations = 10000;
    params.opt.fx_tolerance = 0; // set to 0 to avoid expensive evaluations
    params.opt.dx_tolerance = 1e-9;
    params.opt.verbosity = 0; // silent execution     
    params.opt.use_ell_infinity = 0; 
    params.opt.layer_upper_bound = 0;
    params.opt.use_temporal_penalty = 0;
    params.opt.fx_truth = 0;
    params.opt.fx_truth_eps = 0;    
    params.opt.use_both_penalties = 0;
    params.prob.solve_pixelwise = 1;
    params.opt.use_initial_layers = 0;
    return params;
}

void print_primaldual_params(const primaldual_params& params) {
    mexPrintf("sigma y: %f\n", params.opt.sigma_y);
    mexPrintf("sigma c: %f\n", params.opt.sigma_c);
    mexPrintf("nodes: %d, edges: %d\n", params.prob.nnodes, params.prob.nedges);
    mexPrintf("rows: %d, cols: %d\n", params.prob.rows, params.prob.cols);
    
    mexPrintf("use_ell_infinity : %d\n", params.opt.use_ell_infinity);
    mexPrintf("use_both_penalties : %d\n", params.opt.use_both_penalties);    
    mexPrintf("layer_upper_bound: %f\n", params.opt.layer_upper_bound);
    mexPrintf("use_temporal_penalty (fg prior): %d\n", params.opt.use_temporal_penalty);
    float mu_tau1 = accumulate( params.prob.tau1.begin(), 
                    params.prob.tau1.end(), 0.0f) / params.prob.tau1.size();
    mexPrintf("tau1 (L1): %f tau2 (L-inf): %f\n", mu_tau1, params.prob.tau2);
    mexPrintf("weights size: %d\n", params.prob.weights.size() );
    mexPrintf("occweights size: %d\n", params.prob.occweights.size() );
    mexPrintf("solve pixelwise: %d\n", params.prob.solve_pixelwise );
    mexPrintf("use initialization: %d\n", params.opt.use_initial_layers );
    
    mexPrintf("\n");
    mexPrintf("tau1: %d\n", params.prob.tau1.size() );
    mexPrintf("occweights: %d\n", params.prob.occweights.size() );
    mexPrintf("Wx: %d\n", params.prob.weights.size() );
    mexPrintf("kappa: %d\n", params.prob.kappa.size() );
    mexPrintf("kappa: %d\n", params.prob.kappa.size() );
}