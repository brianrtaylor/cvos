#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std;

int inline linear_index(int r, int c, int k, int rows, int cols);

// USAGE: write_evaluation_track_file_mex( segmentation, out_fname );
// Saves segmentation into an ASCII file, with same format as what is
// used in the MOSEG evaluation software.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // ----------------------------------------------------------------------    
  double* S   = (double*)mxGetPr( prhs[0] );
  
  mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
  int ndims = mxGetNumberOfDimensions( prhs[0] );
  if (ndims < 2) {
    mexPrintf("First argument size is unexpected. Quitting.\n");
    return;
  }
  int rows = sz[0];
  int cols = sz[1];
  int T = sz[2];

  int m = mxGetNumberOfElements( prhs[1] )+1;
  char fname[512];
  if (m > 512) { mexPrintf("crap..\n"); return; }
  mxGetString(prhs[1], fname, m);
  
  FILE* fid = fopen(fname,"w");
  if (fid == 0) { 
      mexPrintf("file: %s\n", fname);
      mexPrintf("could not open file.\n");
      return;
  }
  mexPrintf("%d %d %d\n", rows, cols, T);
  
  int n = 0;
  for (int i = 0; i < T*rows*cols; ++i) { n += !mxIsNaN(S[i]); }
  
  fprintf(fid, "%d\n%d\n", T, n );
  for (int tt = 0; tt < T; ++tt) { 
      mexPrintf(".");
      for (int rr = 0; rr < rows; ++rr) {
          for (int cc = 0; cc < cols; ++cc) { 
              int q = linear_index(rr,cc,tt,rows,cols);
              double val = S[q];
              if (!mxIsNaN(val)) {
                  fprintf(fid, "%d %d\n", (int)val, 1);
                  fprintf(fid, "%d %d %d\n", cc, rr, tt );
              }
          }
      }
  }
  fclose(fid);
  mexPrintf("\n");
}

int inline linear_index(int r, int c, int k, int rows, int cols) {
    return r + c*rows + k*rows*cols;
}