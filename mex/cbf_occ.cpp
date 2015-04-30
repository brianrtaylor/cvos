//----------------------------------------------------------------------------
// Inputs:
// @param: res: residual
// @param: uvf: forward flow t to t+1
// @param: uvb: backward flow t to t-1
// @param: uvf_rev: reverse forward flow t+1 to t
// @param: img: image grayscale
// @param: W: window size (determines sigma_loc) - cuz I need constant size 
//  sliding window I think for histogram version (max 50, but % determined 
//  on image size)
// @param: sigma_loc: sigma for location (sigma_d) (maybe determined by code)
// @param: sigma_flow: simga for flow (sigma_r)
// @param: sigma_res: sigma for residual (new)
// @param: sigma_img: sigma for image intensity (new)
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
    mexPrintf("cbf_occ.\n");
  // ----------------------------------------------------------------------    
  // values for lookup table for exponential.
  int EXP_NUM_ = 2000;
  double EXP_MAX_ = 5.0f;
  // ----------------------------------------------------------------------  
  // declare variables
  double* res_;
  
  double* res = (double*)mxGetPr( prhs[0] );
  double* uvf = (double*)mxGetPr( prhs[1] );
  double* uvb = (double*)mxGetPr( prhs[2] );
  // double* uvf_rev = (double*)mxGetPr( prhs[3] );
  double* img = (double*)mxGetPr( prhs[4] );

  double w_ = (double)mxGetScalar( prhs[5] );
  double sigma_loc = (double)mxGetScalar( prhs[6] );
  double sigma_flow = (double)mxGetScalar( prhs[7] );
  // double sigma_res = (double)mxGetScalar( prhs[8] );
  double sigma_img = (double)mxGetScalar( prhs[9] );

  mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
  int ndims = mxGetNumberOfDimensions( prhs[0] );
  if (ndims < 2) {
    mexPrintf("First argument size is unexpected. Quitting.\n");
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    res_ = (double*)mxGetPr(plhs[0]);
    return;
  }
  int rows = sz[0];
  int cols = sz[1];

  // plhs[0] = mxCreateNumericArray(ndims, sz, mxDOUBLE_CLASS, mxREAL);
  plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
  res_ = (double*)mxGetPr(plhs[0]); 
  
  double sigma_loc_sq_inv = 1.0/(2.0*sigma_loc*sigma_loc);
  double sigma_flow_sq_inv = 1.0/(2.0*sigma_flow*sigma_flow);
  // double sigma_res_sq_inv = 1.0/(2.0*sigma_res*sigma_res);
  double sigma_img_sq_inv = 1.0/(2.0*sigma_img*sigma_img);

  // mexPrintf("rows: %d, cols: %d, chan: %d\n", rows, cols, chan);
  memcpy(res_, res, sizeof(double)*rows*cols);
  
  // create an exponential function lookup table
  double* exp_ = new double[EXP_NUM_];
  for (int i = 0; i < EXP_NUM_; ++i) { 
    double val = double(i)/double(EXP_NUM_)*EXP_MAX_;
    exp_[i] = exp(-val);
  }

  int w = w_;
  int ww = min((int) ceil(sqrt(EXP_MAX_/sigma_loc_sq_inv)), int(w));      
  double gstep = double(w) / double(ww);
          
  int idx, idx1, idx2, idx3, widx, widx1, widx2, widx3;
  for (int x = 0; x < cols; ++x) { 
    int x_min = max(x-ww, 0);
    int x_max = min(x+ww, cols-1);

    for (int y = 0; y < rows; ++y) {
      idx = linear_index(y, x, rows);
      idx1 = linear_index(y, x, 0, rows, cols);
      idx2 = linear_index(y, x, 1, rows, cols);
      idx3 = linear_index(y, x, 2, rows, cols);

      // double r_o = res[ idx ];
      double fu_o = uvf[ idx1 ];
      double fv_o = uvf[ idx2 ];
      double bu_o = uvb[ idx1 ];
      double bv_o = uvb[ idx2 ];
      // double fru_o = uvf_rev[ idx1 ];
      // double frv_o = uvf_rev[ idx2 ];
      double ir_o = img[ idx1 ];
      double ig_o = img[ idx2 ];
      double ib_o = img[ idx3 ];

      // get local window bounds
      int y_min = max(y-ww, 0);
      int y_max = min(y+ww, rows-1);

      double val_res = 0.0;
      double filter_norm = 1e-6;
      for (int cc=x_min; cc <= x_max; ++cc) { 
        for (int rr=y_min; rr <= y_max; ++rr) { 
          widx = linear_index(rr, cc, rows);
          widx1 = linear_index(rr, cc, 0, rows, cols);
          widx2 = linear_index(rr, cc, 1, rows, cols);
          widx3 = linear_index(rr, cc, 2, rows, cols);

          // double r = res[ widx ];
          double fu = uvf[ widx1 ];
          double fv = uvf[ widx2 ];
          double bu = uvb[ widx1 ];
          double bv = uvb[ widx2 ];
          // double fru = uvf_rev[ widx1 ];
          // double frv = uvf_rev[ widx2 ];
          double ir = img[ widx1 ];
          double ig = img[ widx2 ];
          double ib = img[ widx3 ];

          // weight based on flow magnitude difference
          double uvf_dist = ((fu-fu_o)*(fu-fu_o)+(fv-fv_o)*(fv-fv_o))*sigma_flow_sq_inv;
          double uvb_dist = ((bu-bu_o)*(bu-bu_o)+(bv-bv_o)*(bv-bv_o))*sigma_flow_sq_inv;
          // double uvf_rev_dist = ((fru-fru_o)*(fru-fru_o)+(frv-frv_o)*(frv-frv_o))*sigma_flow_sq_inv;
          // weight based on location distance
          double location_dist = ((cc-x)*(cc-x)+(rr-y)*(rr-y))*sigma_loc_sq_inv;
          // weight based on image distance
          double img_dist = ((ir-ir_o)*(ir-ir_o)+(ig-ig_o)*(ig-ig_o)+(ib-ib_o)*(ib-ib_o))*sigma_img_sq_inv;
          // weight based on residual distance
          // double res_dist = ((r-r_o)*(r-r_o))*sigma_res_sq_inv;

          // so slow:
          //double g = exp(-(flow_dist + geometric_dist));
          double d = uvb_dist + uvf_dist + img_dist + location_dist;
          // double d = uvb_dist + uvf_rev_dist; // + location_dist + img_dist;
          // double d = uvb_dist + location_dist + img_dist + res_dist;
          // double d = uvb_dist + location_dist;
          int exp_idx = min( int( ((d * gstep)/EXP_MAX_) *EXP_NUM_), EXP_NUM_-1);
          double g = exp_[exp_idx];

          // mexPrintf("dist: %f  eid: %d   g: %f\n", d, exp_idx, g);

          filter_norm += g;
          val_res += res[widx]*g;
        }
      }
      val_res /= filter_norm;
      res_[idx1] = val_res;
      
    }
  }
  // fin.
  delete exp_;
  mexPrintf("done.\n");
}