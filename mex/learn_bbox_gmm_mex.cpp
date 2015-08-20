//----------------------------------------------------------------------------
// learn_bbox_gmm_mex
//
// learns a gmm color model for pixels within a local shape classifier window
// 
// @param: I0 (MxNx3): RGB image to learn color models from
// @param: layers (MxNxL): layer segmentation from current frame
// @param: boxes: matlab struct array of local shape classifers
// @param: NUM_GMM_CLUSTERS: number of gaussians in the mixture model
// @param: NUM_GMM_REPETITIONS: number of times to run gmm learning
// @param: NUM_GMM_ITERATIONS: number of iterations in run gmm learning
// @param: FGTHRESH: foreground threshold
// @param: BGTHRESH: background threshold
//
// matlab call:
// [mu_fg, cov_fg, pi_fg, mu_bg, cov_bg, pi_bg] = learn_bbox_gmm_mex( ...
//   I0, layers, boxes, NUM_GMM_CLUSTERS, NUM_GMM_REPETITIONS, ...
//   NUM_GMM_ITERATIONS, FGTHRESH, BGTHRESH)
//----------------------------------------------------------------------------
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include "gmm.h" // vlfeat
#include "cvos_common.h"

#define LOG2PI 1.837877066409345

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  //--------------------------------------------------------------------------
  // process inputs and init outputs
  //--------------------------------------------------------------------------
  double* I0              = (double*) mxGetPr( prhs[0] );
  double* layers          = (double*) mxGetPr( prhs[1] );
  const mxArray* boxes    = prhs[2];
  int NUM_GMM_CLUSTERS    = (int) mxGetScalar( prhs[3] );
  int NUM_GMM_REPETITIONS = (int) mxGetScalar( prhs[4] );
  int NUM_GMM_ITERATIONS  = (int) mxGetScalar( prhs[5] );
  double FGTHRESH         = (double) mxGetScalar( prhs[6] );
  double BGTHRESH         = (double) mxGetScalar( prhs[7] );
    
  mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
  int rows = sz[0];
  int cols = sz[1];
  int chan = sz[2];
  if (chan != 3) { mexPrintf("not a 3-channel image\n"); return; }

  int n_boxes = mxGetM(prhs[2]);

  mwSize out_sz[3];
  out_sz[0] = 3;
  out_sz[1] = NUM_GMM_CLUSTERS;
  out_sz[2] = n_boxes;

  plhs[0] = mxCreateNumericArray(3, out_sz, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(3, out_sz, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(NUM_GMM_CLUSTERS, n_boxes, mxREAL);

  plhs[3] = mxCreateNumericArray(3, out_sz, mxDOUBLE_CLASS, mxREAL);
  plhs[4] = mxCreateNumericArray(3, out_sz, mxDOUBLE_CLASS, mxREAL);
  plhs[5] = mxCreateDoubleMatrix(NUM_GMM_CLUSTERS, n_boxes, mxREAL);

  double* fg_mu  = (double*)mxGetPr(plhs[0]);
  double* fg_cov = (double*)mxGetPr(plhs[1]); 
  double* fg_pi  = (double*)mxGetPr(plhs[2]);   
  double* bg_mu  = (double*)mxGetPr(plhs[3]);
  double* bg_cov = (double*)mxGetPr(plhs[4]); 
  double* bg_pi  = (double*)mxGetPr(plhs[5]);   

  for (int i = 0; i < 3 * NUM_GMM_CLUSTERS * n_boxes; ++i) { 
    fg_mu[i] = mxGetNaN();
    bg_mu[i] = mxGetNaN();
    fg_cov[i] = mxGetNaN();
    bg_cov[i] = mxGetNaN();
  }
  for (int i = 0; i < NUM_GMM_CLUSTERS * n_boxes; ++i) {
    fg_pi[i] = mxGetNaN();
    bg_pi[i] = mxGetNaN();
  }

  double center[2];
  int rad, diam, npxmask;
  double *means, *covariances, *priors;

  //--------------------------------------------------------------------------
  // do the work
  //--------------------------------------------------------------------------
  VlGMM* gmm = vl_gmm_new(VL_TYPE_DOUBLE, 3, NUM_GMM_CLUSTERS);
  vl_gmm_set_max_num_iterations(gmm, NUM_GMM_ITERATIONS);
  vl_gmm_set_covariance_lower_bound(gmm, 1e-3);
  vl_gmm_set_initialization(gmm,VlGMMRand);
      
  for (int i = 0; i < n_boxes; ++i) { 
    double* yp = (double*) mxGetPr( mxGetField(boxes, i, "y"));
    double* xp = (double*) mxGetPr( mxGetField(boxes, i, "x"));
    double* rp = (double*) mxGetPr( mxGetField(boxes, i, "r"));

    int my = (int) (yp[0] + 0.5);
    int mx = (int) (xp[0] + 0.5);
    int r  = (int) (rp[0] + 0.5);
    int diam = 2 * r + 1;

    int cy = (int) MAT2C(my);
    int cx = (int) MAT2C(mx);
   
    int ymin = max(0       , (int) (cy - r) );
    int ymax = min(rows - 1, (int) (cy + r) );
    int xmin = max(0       , (int) (cx - r) );
    int xmax = min(cols - 1, (int) (cx + r) );

    int ysz  = ymax - ymin + 1;
    int xsz  = xmax - xmin + 1;
    int numel_box = ysz * xsz;
    size_t n = numel_box * 3;

    vector<double> x1 = vector<double>(); // background
    vector<double> x2 = vector<double>(); // foreground
    x1.reserve(n); 
    x2.reserve(n);

    double* fg_prob = (double*) mxGetPr( mxGetField(boxes, i, "fg_prob"));
    double* bg_prob = (double*) mxGetPr( mxGetField(boxes, i, "bg_prob"));
    double* learn_mask = (double*) mxGetPr( mxGetField(boxes, i,
      "gmm_learn_colour_mask"));

    for (int cc = xmin; cc <= xmax; ++cc) { 
      for (int rr = ymin; rr <= ymax; ++rr) { 
        int idx = linear_index(rr,cc,rows);

        int brr = rr - cy + r;
        int bcc = cc - cx + r;
        int midx = linear_index(brr, bcc, diam); // looks in a full box

        // occluder
        if ((fg_prob[midx] > FGTHRESH) && (learn_mask[midx] == 1)) {
          for (int u = 0; u < chan; ++u) {
            x2.push_back( I0[ linear_index(rr,cc,u,rows,cols) ] );
          }
        } else if ((bg_prob[midx] > BGTHRESH) && (learn_mask[midx] == 1)) {
          for (int u = 0; u < chan; ++u) {
            x1.push_back( I0[ linear_index(rr,cc,u,rows,cols) ] );
          }
        }
      }
    }

    // safe check
    int too_few1 = (x1.size()/chan <= NUM_GMM_CLUSTERS);
    int too_few2 = (x2.size()/chan <= NUM_GMM_CLUSTERS);
    if (too_few1 || too_few2) { continue; }

    //------------------------------------------------------------------------
    // compute gmms
    //------------------------------------------------------------------------
    vl_gmm_cluster(gmm, &x1[0], x1.size()/chan );

    int ind = i * (NUM_GMM_CLUSTERS);
    int ind3 = ind * 3;
    int sz = sizeof(double) * (NUM_GMM_CLUSTERS);
    int sz3 = sz * 3;
    
    means = (double*)vl_gmm_get_means(gmm);
    covariances = (double*)vl_gmm_get_covariances(gmm);
    priors = (double*)vl_gmm_get_priors(gmm);
    
    memcpy(bg_mu  + ind3, means, sz3 );
    memcpy(bg_cov + ind3, covariances, sz3 );
    memcpy(bg_pi  + ind , priors, sz );

    vl_gmm_reset( gmm );
    vl_gmm_cluster(gmm, &x2[0], x2.size()/chan );

    memcpy(fg_mu  + ind3, means, sz3 );
    memcpy(fg_cov + ind3, covariances, sz3 );
    memcpy(fg_pi  + ind , priors, sz );
  }
  vl_gmm_delete(gmm);
}
