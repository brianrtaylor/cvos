//! @file bcv_diff_ops.h
//! A file with code for pixelwise image differences. These
//! operators are used by the TV-optimization algorithms, and some
//! efforts have been made in accelerating them.
//!
//! Forward differences are used throughout ( y(i) = x(i+1)-x(i) ).
//! 'dx' and 'dy' mean horizontal/vertical differences.
//! the suffix 't' implies transpose (i.e. if the difference was
//! represented as a matrix).
//!
//! Circular boundary conditions are assumed.
//!
//! Since MATLAB lays out image differently from OpenCV, a separate
//! set of functions is provided for use with MATLAB. This is done
//! in case, you want to call one of the solvers from MATLAB. 
//! 
//! There are 'designated' functions for 'single-channel' images, and
//! other 'designated' functions for 'multi-channel' images. This is,
//! again, for speed.
#ifndef BCV_DIFF_OPS_H_
#define BCV_DIFF_OPS_H_

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include "bcv_utils.h"

using namespace std;

//! Applies pixelwise gradient operator:
//!         out = [D_x( img ); D_y( img )];
void apply_pixelwise_gradient_op(vector<float>& out, const vector<float>& in, int rows, int cols, int chan);
void apply_pixelwise_gradient_op(vector<float>& out, const vector<float>& in, int rows, int cols);
//! Applies the transpose of the pixelwise gradient operator:
//!         out = D_x^T( edg ) + D_y^T( edg )
void apply_pixelwise_gradient_op_transpose(vector<float>& out, const vector<float>& in, int rows, int cols, int chan); 
void apply_pixelwise_gradient_op_transpose(vector<float>& out, const vector<float>& in, int rows, int cols); 

void apply_dx(float* out, const float* in, int rows, int cols, int chan);
void apply_dxt(float* out, const float* in, int rows, int cols, int chan);
void apply_dy(float* out, const float* in, int rows, int cols, int chan);
//! result is 'ADDED' to the 'out' vector (not 'REPLACED')  
void apply_dyt(float* out, const float* in, int rows, int cols, int chan);

void apply_dx(float* out, const float* in, int rows, int cols);
void apply_dxt(float* out, const float* in, int rows, int cols);
void apply_dy(float* out, const float* in, int rows, int cols);
//! result is 'ADDED' to the 'out' vector (not 'REPLACED')  
void apply_dyt(float* out, const float* in, int rows, int cols);

#if defined(HAVE_SSE) && !defined(HAVE_MATLAB)
void apply_dx_sse(float* out, const float* in, int rows, int cols);
void apply_dxt_sse(float* out, const float* in, int rows, int cols);
void apply_dy_sse(float* out, const float* in, int rows, int cols);
//! result is 'ADDED' to the 'out' vector (not 'REPLACED')    
void apply_dyt_sse(float* out, const float* in, int rows, int cols);
#endif

// because of MATLAB-esque memory layout, it is not necessary to have
// 'single-channel' and 'multi-channel' versions.
#if  !defined(HAVE_SSE) && defined(HAVE_MATLAB)
void apply_dx_matlab(float* out, const float* in, int rows, int cols, int chan=1);
void apply_dxt_matlab(float* out, const float* in, int rows, int cols, int chan=1);
void apply_dy_matlab(float* out, const float* in, int rows, int cols, int chan=1);
//! result is 'ADDED' to the 'out' vector (not 'REPLACED')  
void apply_dyt_matlab(float* out, const float* in, int rows, int cols, int chan=1);
#endif
#if defined(HAVE_SSE) && defined(HAVE_MATLAB)
void apply_dx_matlab_sse(float* out, const float* in, int rows, int cols, int chan=1);
void apply_dxt_matlab_sse(float* out, const float* in, int rows, int cols, int chan=1);
void apply_dy_matlab_sse(float* out, const float* in, int rows, int cols, int chan=1);
//! result is 'ADDED' to the 'out' vector (not 'REPLACED')  
void apply_dyt_matlab_sse(float* out, const float* in, int rows, int cols, int chan=1);
#endif

#endif // SPARSE_OP_H_
