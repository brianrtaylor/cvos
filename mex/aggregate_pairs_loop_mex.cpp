#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std;


#define MAT2C(x) ((x)-1)
#define C2MAT(x) ((x)+1)

void print_usage() { 
    mexPrintf("[unique_pairs, score_pairs, ind_pairs] = ...\n");
    mexPrintf("   aggregate_pairs_loop(pairs, unique_pairs, pair_weights, sorti )\n");
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // ------------------------------------------------------------------------
    // process inputs and init outputs
    // ------------------------------------------------------------------------
    // error checking.
    if (nrhs != 4) { mexPrintf("error.\n"); print_usage(); return; }
    if (nlhs != 3) { mexPrintf("error.\n"); print_usage(); return;  }
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || 
        !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3])) { 
        mexPrintf("arguments should be doubles.\n");
        return; 
    }
    // get input args.
    double* pairs = (double*) mxGetPr( prhs[0] );
    double* unique_pairs = (double*) mxGetPr( prhs[1] );
    double* pair_weights = (double*) mxGetPr( prhs[2] );
    double* sorti = (double*) mxGetPr( prhs[3] );

    int npairs = (int)mxGetM( prhs[0] );
    int nunique_pairs = (int)mxGetM( prhs[1] );

    // ------------------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(nunique_pairs, 2, mxREAL);
    double* out_unique_pairs = (double*)mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(nunique_pairs, 1, mxREAL);
    double* out_score_pairs = (double*)mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(npairs, 1, mxREAL);
    double* out_ind_pairs = (double*)mxGetPr(plhs[2]);
    // ------------------------------------------------------------------------
    memcpy(out_unique_pairs, unique_pairs, sizeof(double)*2*nunique_pairs);
    memset(out_score_pairs, 0, sizeof(double)*nunique_pairs);
    memset(out_ind_pairs, 0, sizeof(double)*nunique_pairs);
    //
    int k = 0; int sk = 0;
    while (k < npairs) {
        // if pairs and unique_pairs match, add pair score to compressed unique_pairs
        // and increment to see next pair. Else, increment to next unique_pair
        if ((pairs[k] == unique_pairs[sk]) && 
            (pairs[k+npairs] == unique_pairs[sk+nunique_pairs]) ) {

            out_score_pairs[sk] = max( out_score_pairs[sk], pair_weights[k] );            
            out_ind_pairs[ MAT2C((int)sorti[k]) ] = C2MAT(sk);
            k+=1;
        } else { sk+=1; }
    }
}