#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <limits> 
#include "cvos_common.h"
// USAGE: mag = occweight_edge_crossing(img, constraints)
//
//  img - edge image; "non-edges" are 'high' and "edges" are LOW (~0)
//  constraints - Nx2 array of indices into image.
//  for each constraint, this function finds the smallest value of 'img'
//  along the line connecting two points.
//
using namespace std;

vector<int> bresenham(int x0, int y0, int x1, int y1, int rows);

#define MAT2C(x) ((x)-1)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if ((nrhs != 2) || (nlhs != 1)) {
    mexPrintf("Incorrect number of input/output arguments.\n");
    mexPrintf("USAGE: mag = occweight_edge_crossing(img, constraints)\n");    
    return;
  }
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) { 
    mexPrintf("MAKE SURE that both arguments are double (they are not currently.\n");
    return;
  }
  // ----------------------------------------------------------------------    
  double* img   = (double*)mxGetPr( prhs[0] );
  double* constraints = (double*)mxGetPr( prhs[1] );

  mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
  int ndims = mxGetNumberOfDimensions( prhs[0] );
  int rows = sz[0];
  int cols = sz[1];
  if ((ndims>2) && (sz[2]>1)) { 
      mexPrintf("expected a single-channel edge map as first argument.\n");
      return;
  }
  int ncons = mxGetM( prhs[1] );
  int two = mxGetN( prhs[1] );
  if (two != 2) {
      mexPrintf(" expected a Nx2 array as the second argument (constraints)\n");
  }
  // ----------------------------------------------------------------------
  // lift-off...
  plhs[0] = mxCreateDoubleMatrix(ncons, 1, mxREAL);
  double* out = (double*)mxGetPr(plhs[0]);
  
  for (int i = 0; i < ncons; ++i) { 
    int i1 = MAT2C( constraints[i] );
    int i2 = MAT2C( constraints[i+ncons] );

    int r1 = getrow(i1, rows);
    int c1 = getcol(i1, rows);
    int r2 = getrow(i2, rows);
    int c2 = getcol(i2, rows);

    vector<int> idx = bresenham(c1, r1, c2, r2, rows);

    out[i] = numeric_limits<double>::max();
    for (size_t j = 0; j < idx.size(); ++j) { 
        out[i] = min(out[i], img[ idx[j] ]);
    }
    //mexPrintf("(%d,%d) (%d,%d) %d - %f\n", c1, r1, c2, r2, idx.size() , out[i] );    
  }
}

// this is from rosetta-code
// (http://rosettacode.org/wiki/Bitmap/Bresenham's_line_algorithm)
vector<int> bresenham(int x0, int y0, int x1, int y1, int rows) {
  vector<int> res = vector<int>();
  res.reserve(200);
    
  int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
  int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1; 
  int err = (dx>dy ? dx : -dy)/2, e2;
 
  for(;;){
    res.push_back( linear_index(y0,x0, rows) );
    if (x0==x1 && y0==y1) break;
    e2 = err;
    if (e2 >-dx) { err -= dy; x0 += sx; }
    if (e2 < dy) { err += dx; y0 += sy; }
  }
  return res;
}