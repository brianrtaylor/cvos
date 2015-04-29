#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <cassert>

using namespace std;
int inline linear_index(int r, int c, int k, int rows, int cols);
int inline linear_index(int r, int c, int rows);
double normpdf(double* I, double* mu, double* cov, int dim);
double utils_eval_gmm(double* I, double* mu, double* cov, double* pi, int dim, int K);
double convert_to_log_ratio(double p1, double p2);
void force_inside(double* pt, int rows, int cols);
double compute_gmm_weight(double p_occr_fg, double p_occr_bg, double p_occd_fg, double p_occd_bg);
double utils_eval_gmm_fast(double* I, double* mu, double* cov, double* pi, int dim, int K);

// log(2 * pi)
#define LOG2PI 1.837877066409345
#define PI 3.141592653589793

#define MAT2C(x) ((x)-1)
#define C2MAT(x) ((x)+1)

//------------------------------------------------------------------------
// [box_fg, box_bg, box_conf] = eval_bboxes_gmm_mex(im2double(I1), bboxes)
//------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  //----------------------------------------------------------------------
  // process inputs and init outputs
  //----------------------------------------------------------------------
  double* I0 = (double*) mxGetPr( prhs[0] );
  const mxArray* boxes = prhs[1];
 
  // TODO: check if n_boxes is 0
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

  double valid[n_boxes];
  for (int k = 0; k < n_boxes; ++k) {
    valid[k] = 1;
  }

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
    double* fg_prob = (double*) mxGetPr( fg_prob_box);
    double* bg_prob = (double*) mxGetPr( bg_prob_box);
    mxSetField(box_fg, i, "prob", fg_prob_box);
    mxSetField(box_bg, i, "prob", bg_prob_box);

    for (int ii = 0; ii < numel_box; ++ii) { 
      fg_prob[ii] = mxGetNaN();
      bg_prob[ii] = mxGetNaN();
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

    // mexPrintf("box: %d, gmm_dim: %d, K: %d\n", i, gmm_dim, K);

    int bad1 = (mxIsNaN(fg_mu[0]) || mxIsNaN(bg_mu[0]));
    int bad2 = (mxIsNaN(fg_cov[0]) || mxIsNaN(bg_cov[0]));
    int bad3 = (mxIsNaN(fg_pi[0]) || mxIsNaN(bg_pi[0]));
    
    if (bad1 || bad2 || bad3) { valid[i] = 0; continue; }

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

        // double p_occr_fg = utils_eval_gmm_fast( I_spot, 
        //   fg_mu, fg_cov, fg_pi, gmm_dim, K);
        // double p_occr_bg = utils_eval_gmm_fast( I_spot, 
        //   bg_mu, bg_cov, bg_pi, gmm_dim, K); 
        // double p_occd_fg = utils_eval_gmm_fast( I_spot, 
        //   fg_mu, fg_cov, fg_pi, gmm_dim, K); 
        // double p_occd_bg = utils_eval_gmm_fast( I_spot, 
        //   bg_mu, bg_cov, bg_pi, gmm_dim, K); 

        // double p_pxl_fg = utils_eval_gmm_fast( I_spot, 
        //   fg_mu, fg_cov, fg_pi, gmm_dim, K);
        // double p_pxl_bg = utils_eval_gmm_fast( I_spot, 
        //   bg_mu, bg_cov, bg_pi, gmm_dim, K); 
        double p_pxl_fg = utils_eval_gmm_fast( I_spot, 
          fg_mu, fg_cov, fg_pi, gmm_dim, K);
        double p_pxl_bg = utils_eval_gmm_fast( I_spot, 
          bg_mu, bg_cov, bg_pi, gmm_dim, K); 

        double p_fg = p_pxl_fg / (p_pxl_fg + p_pxl_bg);
        double p_bg = p_pxl_bg / (p_pxl_fg + p_pxl_bg);

        // double p_fg = convert_to_log_ratio(p_pxl_fg, p_pxl_bg);
        // double p_bg = convert_to_log_ratio(p_pxl_bg, p_pxl_fg);

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

//------------------------------------------------------------
// Helper functions
//------------------------------------------------------------
double compute_gmm_weight(double p_occr_fg, double p_occr_bg, double p_occd_fg, double p_occd_bg) { 
    double P_occr = p_occr_fg / (p_occr_fg + p_occr_bg);
    double P_occd = p_occd_bg / (p_occd_fg + p_occd_bg);
    return (P_occr+P_occd)/2.0;    
}

void force_inside(double* pt, int rows, int cols) { 
    pt[0] = min(max(0.0 , pt[0]), (double)(cols-1) );
    pt[1] = min(max(0.0 , pt[1]), (double)(rows-1) );
} 

int inline linear_index(int r, int c, int k, int rows, int cols) {
    return r + c*rows + k*rows*cols;
}

int inline linear_index(int r, int c, int rows) {
    return r + c*rows;
}

// double utils_eval_gmm_fast(double* I, double* mu, double* cov, double* pi, int dim, int K) {
//     double p = 0;
//     double sqrtcov, f;
//     double *MU, *COV;
//     for (int k = 0; k < K; ++k) {
//         MU = mu + dim*k;
//         COV = cov + dim*k;     
//         // double logp = -dim*LOG2PI; // think this works if  it was dim/2
//         double logp = -log( pow(2 * PI, (0.5 * (double) dim))); // more correct?
//         for (int i = 0; i < dim; ++i) { 
//             sqrtcov = sqrt(COV[i]);
//             f = (I[i]-MU[i])/sqrtcov;
//             logp += (-0.5*f*f - log(sqrtcov));
//         }
//         p += pi[k]*exp(logp);
//     }
//     return p;
// }

double utils_eval_gmm_fast(double* I, double* mu, double* cov, double* pi, int dim, int K) {
  // double p = 0.0f;
  double p2 = 0.0f;
  // double sqrtcov, f;
  double *MU, *COV;
  for (int k = 0; k < K; ++k) {
    MU = mu + dim*k;
    COV = cov + dim*k;     
    double logp = -dim*LOG2PI; // think this works if  it was dim/2
    // double logp = -log( pow(2 * PI, (0.5 * (double) dim))); // more correct?
    double kp = 1.0;
    for (int i = 0; i < dim; ++i) { 
      kp *= (1.0 / sqrt(COV[i] * 2 * PI));
      kp *= exp((-0.5 * (I[i] - MU[i])*(I[i] - MU[i])) / COV[i]);

      // sqrtcov = sqrt(COV[i]);
      // f = (I[i]-MU[i])/sqrtcov;
      // logp += (-0.5*f*f - log(sqrtcov));
    }
    // p += pi[k]*exp(logp);
    p2 += pi[k]*kp;
  }
  // assert(abs(p - p2) < 1e-8);
  return p2;
}


double utils_eval_gmm(double* I, double* mu, double* cov, double* pi, int dim, int K) {
    double p = 0;
    for (int k = 0; k < K; ++k) { 
        double np = normpdf(I, mu + dim*k, cov + dim*k, dim);
        p += pi[k]*np;
    }
    return p;
}

double normpdf(double* I, double* mu, double* cov, int dim) {
    double logp = -dim*LOG2PI; // think this works if  it was dim/2
    // double logp = -log( pow(2 * PI, (0.5 * (double) dim)));
    double sqrtcov, f;
    for (int i = 0; i < dim; ++i) { 
        sqrtcov = sqrt(cov[i]);
        f = (I[i]-mu[i])/sqrtcov;
        logp += (-0.5*f*f - log(sqrtcov));
    }
    return exp(logp);
}

double convert_to_log_ratio(double p1, double p2) {
    p1 = max(min(300.0, log(p1)), -300.0);
    p2 = max(min(300.0, log(p2)), -300.0);
    return (p1-p2);
}
