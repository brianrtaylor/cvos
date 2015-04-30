//----------------------------------------------------------------------------
// USAGE: region = occ_expand(occ, uv)
// Inputs:
// @param: A: occ
// @param: B: uv 
//----------------------------------------------------------------------------

#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include "cvos_common.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    if ((nrhs !=2) || (nlhs!=1)) { 
        mexPrintf("incorrect number of arguments.\n");
        mexPrintf("USAGE: region = occ_expand(occ, uv) \n");
        return;
    }
    double* occ = (double*)mxGetPr( prhs[0] );
    double* uv = (double*)mxGetPr( prhs[1] );
    double* out;
    
    mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
    int ndims = mxGetNumberOfDimensions( prhs[0] );
    if (ndims < 2) {
        mexPrintf("First argument size is unexpected. Quitting.\n");
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        out = (double*)mxGetPr(plhs[0]);
        return;
    }
    int rows = sz[0];
    int cols = sz[1];

    plhs[0] = mxCreateNumericArray(ndims, sz, mxDOUBLE_CLASS, mxREAL);
    out = (double*)mxGetPr(plhs[0]); 
  
    for (int c = 0; c < cols; ++c) { 
        for (int r = 0; r < rows; ++r) { 
            int id1 = linear_index(r,c,0,rows,cols);
            int id2 = linear_index(r,c,1,rows,cols);
            double uv_mag = sqrt( uv[id1]*uv[id1] + uv[id2]*uv[id2] );
            //
            int Q = min(50, max(2,  (int)(uv_mag*occ[id1]) ));
            double val = 0.0;
            for (int xx = max(0, c-Q); xx<=min(cols-1, c+Q); ++xx) { 
                for (int yy = max(0, r-Q); yy<=min(rows-1, r+Q); ++yy) { 
                    int i = linear_index(yy,xx,rows);
                    val = max(occ[i], val);
                }
            }
            out[id1] = val;
        }
    }
}