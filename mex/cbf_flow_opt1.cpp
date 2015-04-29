//----------------------------------------------------------------------------
// Inputs:
// @param: A: uvf
// @param: B: uvb 
// @param: m_change: pixel-wise weights of how much to change 
//  (currently: 1 - residual)
// @param: residual: residual image / occlusion term (discounts flow from 
//  occluded pixels to be counted into local window)
// @param: W: window size (determines sigma_d) - cuz I need constant size 
//  sliding window I think for histogram version (max 50, but % determined 
//  on image size)
// @param: sigma_r: think I used to determine this based on range of flow 
//  within a window 
//----------------------------------------------------------------------------

#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <sys/time.h>

using namespace std;
// Aout = bflt_mex(A, B, mask_change_, mask_ignore_, w_, sigma_d, sigma_r)

int inline linear_index(int r, int c, int k, int rows, int cols);
int inline linear_index(int r, int c, int rows);

unsigned long diff_us(timeval t1, timeval t2) {
  return ((( t1.tv_sec - t2.tv_sec) * 1000000000) + 
          ( t1.tv_usec - t2.tv_usec));
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // ----------------------------------------------------------------------    
  // double MAXG = 100.0;
  // double MINW = 20.0;
  // // double opts_w = 2.0;
  // double opts_TOOMUCHMASKTHRESH = 0.95;
  // values for lookup table for exponential.
  int EXP_NUM_ = 2000;
  double EXP_MAX_ = 5.0f;
  // ----------------------------------------------------------------------  
  // declare variables
  double* uvf_;
  
  double* uvf = (double*)mxGetPr( prhs[0] );
  double* uvb = (double*)mxGetPr( prhs[1] );
  double* m_change = (double*)mxGetPr( prhs[2] );
  double* resf = (double*)mxGetPr( prhs[3] );
  // double* resb = (double*)mxGetPr( prhs[4] );
   
  double w_ = (double)mxGetScalar( prhs[5] );
  // double sigma_uvf = (double)mxGetScalar( prhs[6] );
  double sigma_uvb = (double)mxGetScalar( prhs[7] );
  double sigma_loc = (double)mxGetScalar( prhs[8] );

  struct timeval t_start, t_end;
  unsigned long diff;
  //----------------------------------------------------------------------
  // pass in exponential function lookup table
  //----------------------------------------------------------------------
  gettimeofday(&t_start, NULL);
  double* exp_mat_ = (double*)mxGetPr( prhs[9] );
  gettimeofday(&t_end, NULL);
  diff = diff_us(t_end, t_start);
  mexPrintf("cbf_flow_opt1: reading exponential %0.6f\n", diff);
  //----------------------------------------------------------------------

  mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
  int ndims = mxGetNumberOfDimensions( prhs[0] );
  if (ndims < 2) {
    mexPrintf("First argument size is unexpected. Quitting.\n");
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    uvf_ = (double*)mxGetPr(plhs[0]);
    return;
  }
  int rows = sz[0];
  int cols = sz[1];
  int chan = sz[2];

  // T: just assume 2 channel matrix is passed for now.
  plhs[0] = mxCreateNumericArray(ndims, sz, mxDOUBLE_CLASS, mxREAL);
  uvf_ = (double*)mxGetPr(plhs[0]); 
  
  // double sigma_uvf_sq_inv = 1.0/(2.0*sigma_uvf*sigma_uvf);
  double sigma_uvb_sq_inv = 1.0/(2.0*sigma_uvb*sigma_uvb);
  // double sigma_loc_sq_inv = 1.0/(2.0*sigma_loc*sigma_loc);
  // double sigma_img_sq_inv = 1.0/(2.0*sigma_img*sigma_img);
  
  memcpy(uvf_, uvf, sizeof(double)*rows*cols*chan);
 
  //----------------------------------------------------------------------
  // create an exponential function lookup table
  //----------------------------------------------------------------------
  gettimeofday(&t_start, NULL);
  double* exp_ = new double[EXP_NUM_];
  for (int i = 0; i < EXP_NUM_; ++i) { 
    double val = double(i)/double(EXP_NUM_)*EXP_MAX_;
    exp_[i] = exp(-val);
  }
  gettimeofday(&t_end, NULL);
  diff = diff_us(t_end, t_start);
  mexPrintf("cbf_flow_opt1: computing exponential %0.6f\n", diff);
  //----------------------------------------------------------------------

  int w = w_;
  double sigma_loc_square_2 = 2.0*sigma_loc*sigma_loc;
  int ww = min((int) ceil(sqrt(EXP_MAX_*sigma_loc_square_2)), int(w));
  double gstep = double(w) / double(ww);

  double sigma_loc_sq_inv = 1.0/double(sigma_loc_square_2);
      
          
  int idx, idx1, idx2, widx, widx1, widx2;
  for (int x = 0; x < cols; ++x) { 
    for (int y = 0; y < rows; ++y) {
      idx = linear_index(y, x, rows);
      idx1 = linear_index(y, x, 0, rows, cols);
      idx2 = linear_index(y, x, 1, rows, cols);

      double m_o = m_change[ idx ];
      // double br_o = resb[ idx ];
      double fr_o = resf[ idx ];
      double bu_o = uvb[ idx1 ];
      double bv_o = uvb[ idx2 ];
      // double fu_o = uvf[ idx1 ];
      // double fv_o = uvf[ idx2 ];
      // double r0 = residual[ idx ];

      // if ((m_o > 0.1f) && ((fr_o - br_o) > 0.1f)) {
      if ((m_o > 0.1f)) {
        int y_min = max(y-ww, 0);
        int y_max = min(y+ww, rows-1);
        int x_min = max(x-ww, 0);
        int x_max = min(x+ww, cols-1);

        double val_x = 0.0;
        double val_y = 0.0;
        double filter_norm = 1e-6;
        for (int cc=x_min; cc <= x_max; ++cc) { 
          for (int rr=y_min; rr <= y_max; ++rr) { 
            widx = linear_index(rr, cc, rows);
            widx1 = linear_index(rr, cc, 0, rows, cols);
            widx2 = linear_index(rr, cc, 1, rows, cols);

            // double m = m_change[ widx ];
            // double br = resb[ widx ];
            double fr = resf[ widx ];
            double bu = uvb [ widx1 ];
            double bv = uvb [ widx2 ];
            // double fu = uvf [ widx1 ];
            // double fv = uvf [ widx2 ];

            // weight based on flow magnitude difference
            double uvb_dist_raw = ((bu-bu_o)*(bu-bu_o)+(bv-bv_o)*(bv-bv_o));
            // double uvf_dist_raw = ((fu-fu_o)*(fu-fu_o)+(fv-fv_o)*(fv-fv_o));
            double uvb_dist = uvb_dist_raw*sigma_uvb_sq_inv;
            // double uvf_dist = uvf_dist_raw*sigma_uvf_sq_inv;

            // weight based on geometric distance
            double loc_dist = ((cc-x)*(cc-x)+(rr-y)*(rr-y))*sigma_loc_sq_inv;
            // so slow:
            //double g = exp(-(flow_dist + geometric_dist));
            // double d = uvb_dist + uvf_dist + loc_dist;
            double d = uvb_dist + loc_dist; // 20140807 works
            int exp_idx = min( int( ((d * gstep)/EXP_MAX_) *EXP_NUM_), EXP_NUM_-1);
            double g = exp_[exp_idx] * (1.0 - fr);

            // mexPrintf("dist: %f  eid: %d   g: %f\n", d, exp_idx, g);

            filter_norm += g;
            val_x += uvf[widx1]*g;
            val_y += uvf[widx2]*g;
          }
        }
        val_x /= filter_norm;
        val_y /= filter_norm;
        uvf_[idx1] = val_x;
        uvf_[idx2] = val_y;
      }
    }
  }
  // fin.
  delete exp_;
}

int inline linear_index(int r, int c, int k, int rows, int cols) {
  return r + c*rows + k*rows*cols;
}

int inline linear_index(int r, int c, int rows) {
  return r + c*rows;
}
