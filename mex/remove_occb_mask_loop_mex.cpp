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

#define MAT2C(x) ((x)-1)
#define C2MAT(x) ((x)+1)

void print_usage() { 
  mexPrintf("[out_msfm, occ_cc_max] = remove_occb_mask_loop_mex( ...\n");
  mexPrintf("  in_msfm, occ_msfm, occ_cc_labels)\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // ------------------------------------------------------------------------
  // process inputs and init outputs
  // ------------------------------------------------------------------------
  // error checking.
  if (nrhs != 3) { mexPrintf("error.\n"); print_usage(); return; }
  if (nlhs > 2) { mexPrintf("error.\n"); print_usage(); return; }
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2])) {
    mexPrintf("arguments should be doubles.\n");
  }
  // get input args.
  double* in_msfm = (double*) mxGetPr( prhs[0] );
  double* occ_msfm = (double*) mxGetPr( prhs[1] );
  double* occ_cc_labels = (double*) mxGetPr( prhs[2] );

  int rows = (int)mxGetM( prhs[0] );
  int cols = (int)mxGetN( prhs[0] );
  int npx = rows * cols;

  // figure out how many components there are (max of occ_msfm)
  int ncc = 0;
  for (int i = 0; i < npx; ++i) {
    ncc = max((int) occ_cc_labels[i], ncc);
  }
  // ------------------------------------------------------------------------
  plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
  double* out_msfm = (double*) mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(ncc + 1 + 1, 1, mxREAL);
  double* occ_cc_max = (double*) mxGetPr(plhs[1]);
  // ------------------------------------------------------------------------

  memcpy(out_msfm, in_msfm, sizeof(double) * npx);
  for (int k = 0; k < ncc + 1; ++k) {
    occ_cc_max[k] = 0;
  }

  // per component count the max occ_msfm score
  for (int i = 0; i < npx; ++i) {
    int l = (int) occ_cc_labels[i];
    occ_cc_max[l] = max(occ_cc_max[l], occ_msfm[i]);
  }
  // do subtraction
  for (int i = 0; i < npx; ++i) {
    int l = (int) occ_cc_labels[i];
    if (l == 0) { continue; }
    out_msfm[i] = in_msfm[i] - (1.0 + 2.0 * occ_cc_max[l]);
  }
}

