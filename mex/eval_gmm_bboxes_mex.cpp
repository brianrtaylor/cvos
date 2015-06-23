#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <cassert>
#include "cvos_common.h"
#include "gmm_utils.h"

using namespace std;

//------------------------------------------------------------------------
// [box_fg, box_bg, box_conf] = eval_bboxes_gmm_mex(im2double(I1), bboxes)
//------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  //----------------------------------------------------------------------
  // process inputs and init outputs
  //----------------------------------------------------------------------
  double* I0 = (double*) mxGetPr( prhs[0] );
  const mxArray* boxes = prhs[1];
  int n_boxes = mxGetM(prhs[1]);

  const char* fields[] = {"prob"};
  plhs[0] = mxCreateStructMatrix(n_boxes, 1, 1, fields); 
  plhs[1] = mxCreateStructMatrix(n_boxes, 1, 1, fields);

  mxArray* box_fg = plhs[0];
  mxArray* box_bg = plhs[1];

  mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
  int rows = sz[0];
  int cols = sz[1];
  int chan = sz[2];
  if (chan != 3) { mexPrintf("not a 3-channel image\n"); return; }

  double I_spot[chan];
  int box_idx;

  //----------------------------------------------------------------------
  // do some work
  //----------------------------------------------------------------------
  for (int i = 0; i < n_boxes; ++i) {
    double* yp = (double*) mxGetPr( mxGetField(boxes, i, "y"));
    double* xp = (double*) mxGetPr( mxGetField(boxes, i, "x"));
    double* rp = (double*) mxGetPr( mxGetField(boxes, i, "r"));

    int y = (int) (yp[0] + 0.5);
    int x = (int) (xp[0] + 0.5);
    int r = (int) (rp[0] + 0.5);
    int diam = 2 * r + 1;
    int numel_box = diam * diam;

    // create output arrays
    mxArray* fg_prob_box = mxCreateDoubleMatrix(diam, diam, mxREAL);
    mxArray* bg_prob_box = mxCreateDoubleMatrix(diam, diam, mxREAL);
    double* fg_prob = (double*) mxGetPr( fg_prob_box );
    double* bg_prob = (double*) mxGetPr( bg_prob_box );
    mxSetField(box_fg, i, "prob", fg_prob_box);
    mxSetField(box_bg, i, "prob", bg_prob_box);

    for (int ii = 0; ii < numel_box; ++ii) { 
      fg_prob[ii] = 0.5; //mxGetNaN();
      bg_prob[ii] = 0.5; //mxGetNaN();
    }

    // get gmm info
    double* bg_mu  = (double*) mxGetPr( mxGetField(boxes, i, "bg_gmm_mu") );
    double* bg_cov = (double*) mxGetPr( mxGetField(boxes, i, "bg_gmm_cov") );
    double* bg_pi  = (double*) mxGetPr( mxGetField(boxes, i, "bg_gmm_pi") );
    
    double* fg_mu  = (double*) mxGetPr( mxGetField(boxes, i, "fg_gmm_mu") );
    double* fg_cov = (double*) mxGetPr( mxGetField(boxes, i, "fg_gmm_cov") );
    double* fg_pi  = (double*) mxGetPr( mxGetField(boxes, i, "fg_gmm_pi") );

    mwSize* sz = (mwSize*) mxGetDimensions( mxGetField(boxes, i, "fg_gmm_mu") );
    int gmm_dim = sz[0];
    int K = sz[1];

    int bad1 = (mxIsNaN(fg_mu[0]) || mxIsNaN(bg_mu[0]));
    int bad2 = (mxIsNaN(fg_cov[0]) || mxIsNaN(bg_cov[0]));
    int bad3 = (mxIsNaN(fg_pi[0]) || mxIsNaN(bg_pi[0]));
    
    if (bad1 || bad2 || bad3) { continue; }

    int cy = (int) MAT2C(y);
    int cx = (int) MAT2C(x);

    int ymin = max(0       , (int) (cy - r) );
    int ymax = min(rows - 1, (int) (cy + r) );
    int xmin = max(0       , (int) (cx - r) );
    int xmax = min(cols - 1, (int) (cx + r) );

    for (int cc = xmin; cc <= xmax; ++cc) {
      for (int rr = ymin; rr <= ymax; ++rr) {
        int idx = linear_index(rr,cc,rows);

        for (int u = 0; u < 3; ++u) { 
          I_spot[u] = I0[ linear_index(rr, cc, u, rows, cols) ];
        }

        double p_pxl_fg = utils_eval_gmm_fast( I_spot, 
          fg_mu, fg_cov, fg_pi, gmm_dim, K);
        double p_pxl_bg = utils_eval_gmm_fast( I_spot, 
          bg_mu, bg_cov, bg_pi, gmm_dim, K); 

        double p_fg = p_pxl_fg / (p_pxl_fg + p_pxl_bg);
        double p_bg = p_pxl_bg / (p_pxl_fg + p_pxl_bg);

        if (mxIsNaN(p_fg) || mxIsNaN(p_bg)) {
            mexPrintf("xyi: %f %f\n", p_pxl_fg, p_pxl_bg);
        }
        
        // write out the result
        int brr = rr - cy + r;
        int bcc = cc - cx + r;
        box_idx = linear_index(brr, bcc, diam);
        fg_prob[box_idx] = p_fg;
        bg_prob[box_idx] = p_bg;
      }
    }
  }
}