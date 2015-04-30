#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include "cvos_common.h"

using namespace std;

struct constraint { int occr; int occd; };
struct cluster { 
    vector<constraint> idx; 
    vector<float> mu;
    vector<int> neighbors; // neighboring clusters
    int valid;
};

vector<int> init_label_mask(int rows, int cols, int group_sz);

//! removes clusters without any points. DOES NOT CHANGE 'neighbors' vector.
void prune_constraint_groups(vector<cluster>& in);

//! update the cluster mean
void update_step(vector<cluster>& in, int rows);

//! assign constraints to nearest clusters
void assign_step(vector<cluster>& groups, vector<int>& assignment,
                                const vector<constraint>& cons, int rows);

//! L2 distance between vectors
float compute_ell2(const vector<float>& x1, const vector<float>& x2);

//! 4x1 representation of a constraint (x1,y1,x2,y2)
vector<float> get_constraint_vec(const constraint& c, int rows);

//! convert double array to a reasonable representation
vector<constraint> init_constraints(double* cons, int m);

//! for every group, store indices of neighboring groups.
void init_cluster_neighbors(vector<cluster>& groups, int rows, int group_sz);

void get_initial_assignment(vector<cluster>& groups, vector<int>& assignment,
                                const vector<constraint>& cons, int rows);
                                
// groups = utils_group_constraints(constraints, imsize, group_size)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    if ((nrhs !=3) || (nlhs!=1)) { 
        mexPrintf("incorrect number of arguments.\n");
        mexPrintf("USAGE: groups = utils_group_constraints(cons, imsize, group_sz)\n");
        return;
    }
    double* cons = (double*)mxGetPr( prhs[0] );
    if (mxGetN(prhs[0])!=2) {
        mexPrintf("constraints should form a Nx2 matrix.\n"); return;
    }
    int num_constraints = (int)mxGetM(prhs[0]);
    vector<constraint> constraints = init_constraints(cons, num_constraints);

    double* imsize = (double*)mxGetPr( prhs[1] );
    int group_sz = mxGetScalar( prhs[2] );
    int rows = imsize[0];
    int cols = imsize[1];
    // now can start to do clustering.
    int max_iters = 0;    
    // ------------------------------------------------------------------//  
    vector<int> labels = init_label_mask(rows, cols, group_sz);
    int num_groups = 1+*max_element(labels.begin(), labels.end());
    
    vector<cluster> groups = vector<cluster>(num_groups);
    vector<int> assignment = vector<int>(num_constraints);
    for (int i = 0; i < num_constraints; ++i) { 
        constraint c = constraints[i];
        int id = labels[c.occr];
        groups[ id ].idx.push_back( c );
        groups[ id ].mu = vector<float>(4, 0);
        groups[ id ].valid = 1;
    }
    prune_constraint_groups(groups); // remove empty groups
    update_step(groups, rows); // compute group center
    get_initial_assignment(groups, assignment, constraints, rows);
    
    plhs[0] = mxCreateDoubleMatrix(assignment.size(), 1, mxREAL);
    double* out = (double*)mxGetPr(plhs[0]); 
    for (int i = 0; i < assignment.size(); ++i) { out[i] = C2MAT( assignment[i] ); }
}

// initializes a latticed image
vector<int> init_label_mask(int rows, int cols, int group_sz) {
    int w = max(2, (int)sqrt(group_sz));
    vector<int> z = vector<int>(rows*cols);
    int label = 0;
    for (int c = 0; c < cols; c+=w) { 
        for (int r = 0; r < rows; r+=w) { 
            for (int c1=c; c1<min(cols,c+w); ++c1) { 
                for (int r1=r;r1<min(rows,r+w);++r1){
                    z[ linear_index(r1,c1,rows) ] = label;
                }
            }
            label++;
        }
    }
    return z;
}

//! init-constraints-mask
void prune_constraint_groups(vector<cluster>& in) {
    vector<cluster> out = vector<cluster>();
    for (int i = 0; i < in.size(); ++i) { 
        if (in[i].idx.size()==0) { continue; }
        out.push_back( in[i] );    
    }
    in = out;
}

//! initialize a vector of constraints given a double matrix (from MATLAB)
vector<constraint> init_constraints(double* cons, int m) {
    vector<constraint> C = vector<constraint>(m);
    for (int i = 0; i < m; ++i) { 
        C[i].occr = MAT2C( cons[i] );
        C[i].occd = MAT2C( cons[i+m] );
    }
    return C;
}

//! update step: mean of each cluster is updated.
void update_step(vector<cluster>& in, int rows) {
    for (int i = 0; i < in.size(); ++i) { 
        int N = in[i].idx.size();
        in[i].valid = (N>0);
        if (N==0) { continue; }
        in[i].mu = vector<float>(4, 0.0f);
        for (int j = 0; j < N; ++j) { 
            vector<float> cvec = get_constraint_vec( in[i].idx[j], rows );
            for (int k = 0; k < 4; ++k) { in[i].mu[k] += cvec[k]; }
        }
        for (int j = 0; j < 4; ++j) { in[i].mu[j]/=N; }
    }
}

vector<float> get_constraint_vec(const constraint& c, int rows) {
    vector<float> v(4);
    v[0] = getcol( c.occr , rows );
    v[1] = getrow( c.occr , rows );
    v[2] = getcol( c.occd , rows );
    v[3] = getrow( c.occd , rows );
    return v;
}

float compute_ell2(const vector<float>& x1, const vector<float>& x2) { 
    return (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) +
           (x1[2]-x2[2])*(x1[2]-x2[2]) + (x1[3]-x2[3])*(x1[3]-x2[3]);
}

//! get initial assignment.
void get_initial_assignment(vector<cluster>& groups, vector<int>& assignment,
                                const vector<constraint>& cons, int rows) {
    // clear current assignments
    for (int i = 0; i < cons.size(); ++i) {
        // for each constraint compute a distance..
        vector<float> cvec = get_constraint_vec( cons[i], rows);

        float min_d = 10000000000;
        int argmin = -1;
        for (int j = 0; j < groups.size(); ++j) {
            float d = compute_ell2( cvec, groups[j].mu );
            if (d < min_d) {
                min_d = d;
                argmin = j;
            }
        }
        if (argmin==-1) { mexPrintf("error.\n"); continue; }
        assignment[i] = argmin;
    }
}