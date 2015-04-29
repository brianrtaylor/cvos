#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>

// a simple bilateral filter
//
// usage img = bilateral_filter_mex(img, sigma_d, sigma_I)

using namespace std;
int inline linear_index(int r, int c, int k, int rows, int cols);
int inline linear_index(int r, int c, int rows);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if ((nrhs != 3) || (nlhs != 1)) {
    mexPrintf("Incorrect number of input/output arguments.\n");
    mexPrintf("USAGE: out = bilateral_filter_mex(img, sigma_d, sigma_I)\n");    
    return;
  }
  // ----------------------------------------------------------------------    
  // values for lookup table for exponential.
  int EXP_NUM_ = 2000;
  double EXP_MAX_ = 8;
  // ----------------------------------------------------------------------  
  // declare variables
  double *img, *out;
  double sigma_d, sigma_I;
  int rows, cols, chan;
  
  img   = (double*)mxGetPr( prhs[0] );
  sigma_d = (double)mxGetScalar( prhs[1] );
  sigma_I = (double)mxGetScalar( prhs[2] );

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
  
  plhs[0] = mxCreateNumericArray(ndims, sz, mxDOUBLE_CLASS, mxREAL);
  out = (double*)mxGetPr(plhs[0]); 
  
  // window half size is 4 sigma (wrt geometric distance)
  // upper bound set arbitrarily to 50, mainly to prevent trolling
  int w = (int)min(50.0, max(3.0, 4*sigma_d) );
  
  // precompute divisions
  double sigma_I_sq_inv = 1.0/(2.0*sigma_I*sigma_I);
  double sigma_d_sq_inv = 1.0/(2.0*sigma_d*sigma_d);
  // precompute exponential
  double* exp_ = new double[EXP_NUM_];
  for (int i = 0; i < EXP_NUM_; ++i) { 
      double val = double(i)/double(EXP_NUM_)*EXP_MAX_;
      exp_[i] = exp(-val);
  }
  //
  int idx1, idx2;
  double* values = new double[chan];
  double* out_values = new double[chan];

  for (int x = 0; x < cols; ++x) {
    // get local window bounds
    int x_min = max(x-w, 0);
    int x_max = min(x+w, cols-1);    
    for (int y = 0; y < rows; ++y) {
        for (int k = 0; k < chan; ++k) { 
            values[k] = img[ linear_index(y, x, k, rows, cols) ];
            out_values[k] = 0.0;
        }
        // get local window bounds
        int y_min = max(y-w, 0);
        int y_max = min(y+w, rows-1);
        
        double filter_norm = 1e-8;
        for (int cc=x_min; cc <= x_max; ++cc) { 
            for (int rr=y_min; rr <= y_max; ++rr) { 
                // weight based on geometric distance
                double geometric_dist = ( (cc-x)*(cc-x) + (rr-y)*(rr-y) )*sigma_d_sq_inv;
                
                // intensity weight
                double dI = 0.0;
                for (int k = 0; k < chan; ++k) { 
                    double val = img[ linear_index(rr, cc, k, rows, cols) ];
                    dI += (val-values[k])*(val-values[k]);
                }
                dI *= sigma_I_sq_inv;
                double val = dI + geometric_dist;
                
                //double g = exp(-(flow_dist + geometric_dist));
                int exp_idx = min( int(val/EXP_MAX_*EXP_NUM_), EXP_NUM_-1);
                if (exp_idx < 0) { exp_idx = EXP_NUM_-1; }
                double g = exp_[exp_idx];
                filter_norm += g;                
                for (int k = 0; k < chan; ++k) { 
                    out_values[k] += g * img[ linear_index(rr, cc, k, rows, cols) ];
                }
            }
        }      
        for (int k = 0; k < chan; ++k) { 
            out_values[k] /= filter_norm;
            out[ linear_index(y, x, k, rows, cols) ] = out_values[k];
        }
    }
  }
  // fin.
  delete exp_;
  delete values;
  delete out_values;
}


int inline linear_index(int r, int c, int k, int rows, int cols) {
    return r + c*rows + k*rows*cols;
}

int inline linear_index(int r, int c, int rows) {
    return r + c*rows;
}