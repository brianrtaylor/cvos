#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>

// Given a flow image, computes the largest difference in angle at each point,
// wrt its neighbors.
// d = flow_angle_distance_mex(uv)

using namespace std;
int inline linear_index(int r, int c, int k, int rows, int cols);
int inline linear_index(int r, int c, int rows);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if ((nrhs != 1) || (nlhs != 1)) {
    mexPrintf("Incorrect number of input/output arguments.\n");
    mexPrintf("USAGE: out = flow_angle_distance_mex(uv)\n");    
    return;
  }
  // ----------------------------------------------------------------------    
  // declare variables
  double *img, *out;
  int rows, cols, chan;
  
  img   = (double*)mxGetPr( prhs[0] );

  mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
  int ndims = mxGetNumberOfDimensions( prhs[0] );
  if (ndims < 2) {
    mexPrintf("First argument size is unexpected. Quitting.\n");
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    out = (double*)mxGetPr(plhs[0]);
    return;
  }
  rows = sz[0];
  cols = sz[1];
  chan = sz[2];
  
  sz[2] = 1;
  plhs[0] = mxCreateNumericArray(ndims, sz, mxDOUBLE_CLASS, mxREAL);
  out = (double*)mxGetPr(plhs[0]); 
  
  int id1, id2;
  double norm, d, d_, theta_x, theta_y, alpha_x, alpha_y;
  for (int x = 0; x < cols; ++x) {
    for (int y = 0; y < rows; ++y) { 
        // get angle:
        id1 = linear_index(y, x, 0, rows, cols);
        id2 = linear_index(y, x, 1, rows, cols);
        norm = sqrt( img[id1]*img[id1] + img[id2]*img[id2] );
        theta_x = img[id1]/norm;
        theta_y = img[id2]/norm;
        d = 0.0;
        for (int xx = max(0, x-1); xx <= min(cols-1,x+1); ++xx) { 
            for (int yy = max(0, y-1); yy <= min(rows-1,y+1); ++yy) { 
                id1 = linear_index(yy, xx, 0, rows, cols);
                id2 = linear_index(yy, xx, 1, rows, cols);
                norm = sqrt( img[id1]*img[id1] + img[id2]*img[id2] );
                alpha_x = img[id1]/norm;
                alpha_y = img[id2]/norm;
                
                d_ = acos( alpha_x*theta_x + alpha_y*theta_y );
                d = max(d, d_);
            }
        }
        id1 = linear_index(y, x, rows);
        out[id1] = d;
    }
  }
}


int inline linear_index(int r, int c, int k, int rows, int cols) {
    return r + c*rows + k*rows*cols;
}

int inline linear_index(int r, int c, int rows) {
    return r + c*rows;
}
