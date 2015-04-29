#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std;
// Aout = bflt_mex(A, B, mask_change_, mask_ignore_, w_, sigma_d, sigma_r)

int inline linear_index(int r, int c, int k, int rows, int cols);
int inline linear_index(int r, int c, int rows);
int determine_region_size(double* A, int r, int c, int w_, 
        int rows, int cols, double MAXG, double MINW, double opts_w);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // ----------------------------------------------------------------------    
  double MAXG = 100.0;
  double MINW = 20.0;
  double opts_w = 2.0;
  double opts_TOOMUCHMASKTHRESH = 0.95;
  // values for lookup table for exponential.
  int EXP_NUM_ = 2000;
  double EXP_MAX_ = 8;
  // ----------------------------------------------------------------------  
  // declare variables
  double *A, *B, *mask_change_, *mask_ignore_, *Aout;
  double w_, sigma_d, sigma_r;
  int rows, cols, chan;
  
  A   = (double*)mxGetPr( prhs[0] );
  B   = (double*)mxGetPr( prhs[1] );
  mask_change_  = (double*)mxGetPr( prhs[2] );
  mask_ignore_  = (double*)mxGetPr( prhs[3] );
   
  w_ = (double)mxGetScalar( prhs[4] );
  sigma_d = (double)mxGetScalar( prhs[5] );
  sigma_r = (double)mxGetScalar( prhs[6] );

  mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
  int ndims = mxGetNumberOfDimensions( prhs[0] );
  if (ndims < 2) {
    mexPrintf("First argument size is unexpected. Quitting.\n");
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    Aout = (double*)mxGetPr(plhs[0]);
    return;
  }
  rows = sz[0];
  cols = sz[1];
  chan = sz[2];
  // T: just assume 2 channel matrix is passed for now.
  plhs[0] = mxCreateNumericArray(ndims, sz, mxDOUBLE_CLASS, mxREAL);
  Aout = (double*)mxGetPr(plhs[0]); 
  
  double sigma_r_sq_inv = 1.0/(2.0*sigma_r*sigma_r);
  double sigma_d_sq_inv = 1.0/(2.0*sigma_d*sigma_d);
  
  memcpy(Aout, A, sizeof(double)*rows*cols*chan);
  
  // create a 
  double* exp_ = new double[EXP_NUM_];
  for (int i = 0; i < EXP_NUM_; ++i) { 
      double val = double(i)/double(EXP_NUM_)*EXP_MAX_;
      exp_[i] = exp(-val);
  }
          
  int mask_idx, idx1, idx2;
  for (int x = 0; x < cols; ++x) { 
    for (int y = 0; y < rows; ++y) {
        mask_idx = linear_index(y, x, rows);
        if (!mask_change_[mask_idx]) { continue; }
        // extract local region from flow magnitude:
        int w = determine_region_size(A, y, x, w_, rows, cols, MAXG, MINW, opts_w);
        // get local window bounds
        int y_min = max(y-w, 0);
        int y_max = min(y+w, rows-1);
        int x_min = max(x-w, 0);
        int x_max = min(x+w, cols-1);
        
        // determine if too masked...
        int num_ignore = 0;
        int num = 0;
        for (int cc = x_min; cc <= x_max; ++cc) { 
            for (int rr=y_min; rr <= y_max; ++rr) { 
                num_ignore += mask_ignore_[ linear_index(rr, cc, rows) ];
                num++;
            }
        }
        int toomasked = num_ignore > (opts_TOOMUCHMASKTHRESH*num);
        if (toomasked) { continue; }        
        
        double bx_o = B[ linear_index(y, x, 0, rows, cols) ];
        double by_o = B[ linear_index(y, x, 1, rows, cols) ];
        
        double val_x = 0.0;
        double val_y = 0.0;
        double filter_norm = 0.0;
        for (int cc=x_min; cc <= x_max; ++cc) { 
            for (int rr=y_min; rr <= y_max; ++rr) { 
                mask_idx = linear_index(rr, cc, rows);
                if (mask_ignore_[mask_idx]) { continue; }
                
                double bx = B[ linear_index(rr, cc, 0, rows, cols) ];
                double by = B[ linear_index(rr, cc, 1, rows, cols) ];
                // weight based on flow magnitude difference
                double flow_dist = ( (bx-bx_o)*(bx-bx_o) + (by-by_o)*(by-by_o) )*sigma_r_sq_inv;
                // weight based on geometric distance
                double geometric_dist = ( (cc-x)*(cc-x) + (rr-y)*(rr-y) )*sigma_d_sq_inv;
                // so slow:
                //double g = exp(-(flow_dist + geometric_dist));
                int exp_idx = min( int( (flow_dist+geometric_dist)/EXP_MAX_*EXP_NUM_), EXP_NUM_-1);
                double g = exp_[exp_idx];
                
                idx1 = linear_index(rr, cc, 0, rows, cols);
                idx2 = linear_index(rr, cc, 1, rows, cols);
                filter_norm += g;
                val_x += A[idx1]*g;
                val_y += A[idx2]*g;
            }
        }
        val_x /= filter_norm;
        val_y /= filter_norm;
        
        idx1 = linear_index(y, x, 0, rows, cols);
        idx2 = linear_index(y, x, 1, rows, cols);
        Aout[idx1] = val_x;
        Aout[idx2] = val_y;
    }
  }
  // fin.
  delete exp_;
}

//! determine region size using maximum flow in a local region
int determine_region_size(double* A, int r, int c, int w_, 
        int rows, int cols, double MAXG, double MINW, double opts_w) {
    int r_min = max(r-w_,0);
    int r_max = min(r+w_, rows-1);
    int c_min = max(c-w_,0);
    int c_max = min(c+w_, cols-1);
    
    double uvx = 0.0;
    double uvy = 0.0;
    
    int n = 0;
    for (int cc = c_min; cc <= c_max; ++cc) { 
        for (int rr = r_min; rr <= r_max; ++rr) {
            int id1 = linear_index(rr, cc, 0, rows, cols);
            int id2 = linear_index(rr, cc, 1, rows, cols);
            double Ax = A[ id1 ];
            double Ay = A[ id2 ];
            if ((Ax==Ax) && (Ay == Ay)) {
                uvx += abs( Ax );
                uvy += abs( Ay );
                n++;
            }
        }
    }
    uvx /= double(n);
    uvy /= double(n);
    
    int maxflow = ceil(max(uvx, uvy));    
    int w = (int)min(MAXG, max(MINW, opts_w*maxflow));
    return w;
}


int inline linear_index(int r, int c, int k, int rows, int cols) {
    return r + c*rows + k*rows*cols;
}

int inline linear_index(int r, int c, int rows) {
    return r + c*rows;
}