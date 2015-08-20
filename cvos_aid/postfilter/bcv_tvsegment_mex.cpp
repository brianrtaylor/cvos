#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>

#include <tvsegment.h>

using namespace std;

// USAGE:
// x = bcv_tvsegment_mex(img, unary, lambda, beta, max_iterations)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // --------------------------------------------------------------------
    if (nrhs !=5) { 
        mexPrintf("Arguments seem to be wrong. Should be:\n");
        mexPrintf("x = bcv_tvsegment_mex(img, unary, lambda, beta, max_iterations)\n");
        return; 
    }
    if (nlhs > 2) { 
        mexPrintf("too many return args.\n");
        return;
    }
    if (!mxIsDouble( prhs[0])) { 
        mexPrintf("argument 1 (image) should be double.\n");
        return;
    }
    if (!mxIsDouble( prhs[1])) { 
        mexPrintf("argument 2 (unary) should be double.\n");
        return;
    }    
    // --------------------------------------------------------------------
    
    
    double* img = (double*)mxGetPr( prhs[0] );
    mwSize* sz = (mwSize*)mxGetDimensions(prhs[0]);
    int ndim = mxGetNumberOfDimensions(prhs[0]);
    int rows = sz[0];
    int cols = sz[1];
    int chan = 1;
    if (ndim > 2) { int chan = sz[2]; }

    double* unary = (double*)mxGetPr( prhs[1] );
    sz = (mwSize*)mxGetDimensions(prhs[1]);
    ndim = mxGetNumberOfDimensions(prhs[1]);
    if (ndim < 3) { mexPrintf("crap"); return; }
    int K = sz[2];

    float lambda = (float)mxGetScalar( prhs[2] );
    float beta = (float)mxGetScalar( prhs[3] );
    int max_iters = (int)mxGetScalar( prhs[4] );

    // same size as the single-channel image
    sz[2] = 1;
    plhs[0] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);
    double* out = (double*)mxGetPr(plhs[0]);
    
    // --------------------------------------------------------------------
    //                      set up the problem
    // --------------------------------------------------------------------
    int n = 1; for (int i = 0; i < ndim; ++i) { n*=sz[i]; }
    tvsegment_params p;
    p.rows = rows;
    p.cols = cols;
    p.chan = chan;
    p.lambda = lambda;
    p.beta = beta;
    p.max_iters = max_iters;
    p.isotropic = 1;
    p.num_clusters = K;
    p.unary = vector<float>(rows*cols*K);
    
    for (size_t i = 0; i < p.unary.size(); ++i) { p.unary[i] = unary[i]; }
    vector<float> imgf = vector<float>(rows*cols*chan);
    for (size_t i = 0; i < imgf.size(); ++i) { imgf[i] = img[i]; }
    // --------------------------------------------------------------------
    //                      solve the problem
    // --------------------------------------------------------------------
    tvsegment tvs = tvsegment(imgf, &p);
    vector<int> seg = tvs.get_assignments();
    
    for (int i = 0; i < rows*cols; ++i) { out[i] = (double)seg[i]; }

    if (nlhs > 1) { 
        vector<float> weights = tvs.get_weights();
        // same size as the single-channel image
        plhs[1] = mxCreateNumericArray(ndim, sz, mxDOUBLE_CLASS, mxREAL);
        double* wout = (double*)mxGetPr(plhs[1]);
        memset(wout, 0, sizeof(double)*rows*cols);
        for (int i = 0; i < rows*cols; ++i) { 
            for (int k = 0; k < K; ++k) { 
                wout[i] += weights[ i + k*rows*cols ];
                wout[i] += weights[ i + (k+K)*rows*cols];                
            }
            wout[i] /= (double(K)*2.0);
        }
    }
}
