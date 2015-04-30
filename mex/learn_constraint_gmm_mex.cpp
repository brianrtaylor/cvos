#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include "gmm.h" // from VLfeat
#include "cvos_common.h"

using namespace std;

#define LOG2PI 1.837877066409345

#define MAT2C(x) ((x)-1)
#define C2MAT(x) ((x)+1)

// [mu_fg, cov_fg, pi_fg, mu_bg, cov_bg, pi_bg] = 
//  learn_constraint_gmm_mex( I0, occd_smooth, occr_smooth, W_binary, PTS1, PTS2, ...
//                    EPS, NUM_GMM_CLUSTERS, NUM_GMM_REPETITIONS, NUM_GMM_ITERATIONS) 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // --------------------------------------------------------------------    
    double* I0   = (double*)mxGetPr( prhs[0] );
    mxLogical* occd_smooth = (mxLogical*)mxGetPr( prhs[1] );
    mxLogical* occr_smooth = (mxLogical*)mxGetPr( prhs[2] );
    mxLogical* W_binary = (mxLogical*)mxGetPr( prhs[3] );
    double* layers = (double*)mxGetPr( prhs[4] );
    double* PTS_OCCD = (double*)mxGetPr( prhs[5] );
    double* PTS_OCCR = (double*)mxGetPr( prhs[6] );
    double* groups = (double*)mxGetPr( prhs[7] );
    
    int EPS = (int)mxGetScalar( prhs[8] );
    int NUM_GMM_CLUSTERS = (int)mxGetScalar( prhs[9] );
    int NUM_GMM_REPETITIONS = (int)mxGetScalar( prhs[10] );
    int NUM_GMM_ITERATIONS = (int)mxGetScalar( prhs[11] );
    
    mwSize* sz = (mwSize*)mxGetDimensions( prhs[0] );
    int rows = sz[0];
    int cols = sz[1];
    int chan = sz[2];
    if (chan != 3) { mexPrintf("not a 3-channel image\n"); return; }

  int ncons = mxGetM(prhs[5]); //max( mxGetM(prhs[5]), mxGetN(prhs[5]) );

  bool LAYER_MASK_GIVEN = (int)max( mxGetM(prhs[4]), mxGetN(prhs[4]) ) > 0;
 
  
  int ngroups = 0;
  for (int i = 0; i < ncons; ++i) { ngroups = max((int)groups[i],ngroups); }
  
  mwSize out_sz[3];
  out_sz[0] = 3;
  out_sz[1] = NUM_GMM_CLUSTERS;
  out_sz[2] = ngroups;
  plhs[0] = mxCreateNumericArray(3, out_sz, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(3, out_sz, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(NUM_GMM_CLUSTERS, ngroups, mxREAL);

  plhs[3] = mxCreateNumericArray(3, out_sz, mxDOUBLE_CLASS, mxREAL);
  plhs[4] = mxCreateNumericArray(3, out_sz, mxDOUBLE_CLASS, mxREAL);
  plhs[5] = mxCreateDoubleMatrix(NUM_GMM_CLUSTERS, ngroups, mxREAL);

  double* fg_mu  = (double*)mxGetPr(plhs[0]);
  double* fg_cov = (double*)mxGetPr(plhs[1]); 
  double* fg_pi  = (double*)mxGetPr(plhs[2]);   
  double* bg_mu  = (double*)mxGetPr(plhs[3]);
  double* bg_cov = (double*)mxGetPr(plhs[4]); 
  double* bg_pi  = (double*)mxGetPr(plhs[5]);   
  // ----------------------------------------------------------------------
  //                                do the work
  // ----------------------------------------------------------------------
  double center[2];
  double *means, *covariances, *priors;
  
  for (int i = 0; i < 3*NUM_GMM_CLUSTERS*ngroups; ++i) { 
    fg_mu[i] = mxGetNaN();
    bg_mu[i] = mxGetNaN();
    fg_cov[i] = mxGetNaN();
    bg_cov[i] = mxGetNaN();
  }
  for (int i = 0; i < NUM_GMM_CLUSTERS*ngroups; ++i) {
    fg_pi[i] = mxGetNaN();
    bg_pi[i] = mxGetNaN();
  }
  

  // ------------------------------------------------------------------
  VlGMM* gmm = vl_gmm_new(VL_TYPE_DOUBLE, 3, NUM_GMM_CLUSTERS);
  vl_gmm_set_max_num_iterations(gmm, NUM_GMM_ITERATIONS);
  vl_gmm_set_covariance_lower_bound(gmm, 1e-3);
  vl_gmm_set_initialization(gmm,VlGMMRand);
      
  for (int grp = 0; grp < ngroups; ++grp) { 
      
      vector<double> x1 = vector<double>();
      vector<double> x2 = vector<double>();
      x1.reserve(3*1024);
      x2.reserve(3*1024);
          
      for (int i = 0; i < ncons; ++i) { 
          if (MAT2C(groups[i]) != grp) { continue; }

          double yy1 = MAT2C( PTS_OCCD[i] );
          double yy2 = MAT2C( PTS_OCCR[i] );
          double xx1 = MAT2C( PTS_OCCD[i+ncons] );
          double xx2 = MAT2C( PTS_OCCR[i+ncons] );
      
          center[0] = (yy1+yy2)/2.0;
          center[1] = (xx1+xx2)/2.0;
          double D = sqrt( (xx2-xx1)*(xx2-xx1)+(yy2-yy1)*(yy2-yy1) );
          //
          int BALL_SIZE = max( EPS, (int)ceil( D/2*1.25) );

          int ball_row1 = max(0, (int)(center[0]-BALL_SIZE) );
          int ball_row2 = min(rows-1, (int)(center[0]+BALL_SIZE) );
          int ball_col1 = max(0, (int)(center[1]-BALL_SIZE) );
          int ball_col2 = min(cols-1, (int)(center[1]+BALL_SIZE) );
      
          if (!LAYER_MASK_GIVEN) { 
          //-------------------------------------------------------------------
              for (int cc = ball_col1; cc <= ball_col2; ++cc) { 
                  for (int rr = ball_row1; rr <= ball_row2; ++rr) { 
                      int idx = linear_index(rr,cc,rows);
                      // put in occluded
                      if (occd_smooth[idx] && !W_binary[idx]) {
                          for (int u = 0; u < chan; ++u) { 
                              x1.push_back( I0[ linear_index(rr,cc,u,rows,cols) ] );
                          }
                      }
                      // put in occluder
                      if (occr_smooth[idx] && !W_binary[idx]) {
                          for (int u = 0; u < chan; ++u) { 
                              x2.push_back( I0[ linear_index(rr,cc,u,rows,cols) ] );
                          }
                      }
                  }
              }
          //-------------------------------------------------------------------
          } else {
          //-------------------------------------------------------------------          
              double occr_layer_value = layers[linear_index(yy2, xx2, rows)];
              double occd_layer_value = layers[linear_index(yy1, xx1, rows)];

              for (int cc = ball_col1; cc <= ball_col2; ++cc) { 
                  for (int rr = ball_row1; rr <= ball_row2; ++rr) { 
                      int idx = linear_index(rr,cc,rows);
                      // put in occluded
                      int OCCR = abs(layers[idx] - occr_layer_value) < abs(layers[idx] - occd_layer_value);

                      if (occd_smooth[idx] && !W_binary[idx] && !OCCR) {
                          for (int u = 0; u < chan; ++u) { 
                              x1.push_back( I0[ linear_index(rr,cc,u,rows,cols) ] );
                          }
                      }
                      // put in occluder
                      if (occr_smooth[idx] && !W_binary[idx] && OCCR) {
                          for (int u = 0; u < chan; ++u) { 
                              x2.push_back( I0[ linear_index(rr,cc,u,rows,cols) ] );
                          }
                      }
                  }
              }
          }
      }
      // ------------------------------------------------------------------
      //                 learn GMM using these points...
      // ------------------------------------------------------------------
      int too_few1 = (x1.size()/chan <= NUM_GMM_CLUSTERS);
      int too_few2 = (x2.size()/chan <= NUM_GMM_CLUSTERS);
      
      if (too_few1 || too_few2) { continue; }
      
      // ------------------------------------------------------------------
      vl_gmm_cluster(gmm, &x1[0], x1.size()/chan );

      means = (double*)vl_gmm_get_means(gmm);
      covariances = (double*)vl_gmm_get_covariances(gmm);
      priors = (double*)vl_gmm_get_priors(gmm);
      
      memcpy(bg_mu  + grp*(NUM_GMM_CLUSTERS*3), means, sizeof(double)*NUM_GMM_CLUSTERS*3 );
      memcpy(bg_cov + grp*(NUM_GMM_CLUSTERS*3), covariances, sizeof(double)*NUM_GMM_CLUSTERS*3 );
      memcpy(bg_pi  + grp*(NUM_GMM_CLUSTERS), priors, sizeof(double)*NUM_GMM_CLUSTERS );
      // ------------------------------------------------------------------

      vl_gmm_reset( gmm );
      vl_gmm_cluster(gmm, &x2[0], x2.size()/chan );

      memcpy(fg_mu  + grp*(NUM_GMM_CLUSTERS*3), means, sizeof(double)*NUM_GMM_CLUSTERS*3 );
      memcpy(fg_cov + grp*(NUM_GMM_CLUSTERS*3), covariances, sizeof(double)*NUM_GMM_CLUSTERS*3 );
      memcpy(fg_pi  + grp*(NUM_GMM_CLUSTERS), priors, sizeof(double)*NUM_GMM_CLUSTERS );
  }
  vl_gmm_delete(gmm);
}