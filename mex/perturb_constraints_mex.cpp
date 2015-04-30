#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include "cvos_common.h"
#include "gmm_utils.h"

using namespace std;

struct gmm { 
    double* mu;
    double* cov;
    double* pi;
    int K;
    int dim;
};

struct constraint_data {
    vector<double> occd_pts;
    vector<double> occr_pts;
    vector<double> weights;
    vector<double> separation;
    gmm gmm_fg;
    gmm gmm_bg;
    int num_pts;
};

struct xform { 
    double scale;
    double theta;
    double T[2];
};

struct memorization_table {
    // stores indices into arrays below:
    vector<short> idx;
    vector<double> data_log;
    vector<double> data_prob;    
    int num;
    int n;
    int rows;
    int cols;
};

struct perturb_params {
  float EDGE_DISTANCE_COST;
  float TRANSLATION_COST;
  int MAX_SHIFT;
  int NUM_SHIFT;
  int NUM_ROTATION;
  float MAX_ROTATION;
  float ROTATION_COST;
  float MAX_STRETCH;
  int NUM_STRETCH;
  float STRETCH_COST;
  float INVALID_COST;    
};

void force_inside(double* pt, int rows, int cols);
void rotate_and_scale( double* out, double* in, double* center, double theta, double scale );

void optimal_constraint_deformation_A2(double* I0, double* edge_image, char* invalid,
        constraint_data* data, vector<int>& shifts_, vector<double>& rotations_, 
        const perturb_params* params, int rows, int cols);

memorization_table init_memorization_table(int rows, int cols, int num_pts);
void inline put_in_memorization_table(memorization_table& m, double log, double prob, int r, int c, int pt);
short inline is_in_memorization_table(memorization_table& m, int r, int c, int pt);
int get_num_groups(const double* groups, int n);
void get_group_assignments(vector<vector<int> >& assignments, double* groups, int ngroups, int n);
int is_gmm_bad(gmm* g);
void adjust_num_stretch(int* NUM_STRETCH, int* MAX_STRETCH, 
                            double separation, double maxallowedseparation);
double get_min_nonzero_separation(double sep, int NUM_STRETCH, int MAX_STRETCH);

// USAGE:    
// [pts_occd_out, pts_occr_out, gmm_weights] = perturb_constraints_mex(...
//               I0, d_edge, pts_occd, pts_occr, pts_max_separation, ...
//               constraint_groups, local_gmm_fg, local_gmm_bg, params );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // ----------------------------------------------------------------------    
  double* I0   = (double*)mxGetPr( prhs[0] );
  double* edge_image   = (double*)mxGetPr( prhs[1] );  
  
  double* pts_occd = (double*)mxGetPr( prhs[2] );
  double* pts_occr = (double*)mxGetPr( prhs[3] );
  double* pts_max_separation = (double*)mxGetPr( prhs[4] );
  double* groups   = (double*)mxGetPr( prhs[5] );
  
  const mxArray* fg_gmm = prhs[6];
  const mxArray* bg_gmm = prhs[7];
  double* gmm_fg_mu  = (double*)mxGetPr( mxGetField(fg_gmm, 0, "mu") );   
  double* gmm_fg_cov = (double*)mxGetPr( mxGetField(fg_gmm, 0, "cov") );
  double* gmm_fg_pi  = (double*)mxGetPr( mxGetField(fg_gmm, 0, "pi") );
    
  double* gmm_bg_mu  = (double*)mxGetPr( mxGetField(bg_gmm, 0, "mu") );
  double* gmm_bg_cov = (double*)mxGetPr( mxGetField(bg_gmm, 0, "cov") );
  double* gmm_bg_pi  = (double*)mxGetPr( mxGetField(bg_gmm, 0, "pi") );
  
  mwSize* sz = (mwSize*) mxGetDimensions( mxGetField(fg_gmm, 0, "mu") );
  int gmm_dim = sz[0];
  int K = sz[1];
  int ncons = mxGetM( prhs[2] ); // get it from constraints input
  
  // various checks to make sure arguments are sane:
  if (ncons != mxGetM(prhs[3])) { 
    mexPrintf("number of occd and occr points doesnt match.\n");
    mexPrintf("check that this behavior is sane.\n");    
  }
  if (ncons != mxGetNumberOfElements(prhs[4])) { 
    mexPrintf("number of constraints doesnt match the number of values in max-separation vector.\n");
    mexPrintf("these numbers should be the same!.\n");
    mexPrintf("check that this behavior is sane.\n");    
  }
  if (ncons != mxGetNumberOfElements(prhs[5])) { 
    mexPrintf("number of constraints doesnt match the number of groups.\n");
    mexPrintf("these numbers should be the same!.\n");
    mexPrintf("check that this behavior is sane.\n");    
  }
  
  
  if (nrhs < 9) { mexPrintf("number of arguments is wrong.\n"); return; }
    
  const mxArray* params_struct = prhs[8];
  perturb_params params;
  params.EDGE_DISTANCE_COST = (float)mxGetScalar(
                        mxGetField(params_struct, 0, "EDGE_DISTANCE_COST"));
  params.TRANSLATION_COST = (float)mxGetScalar(
                        mxGetField(params_struct, 0, "TRANSLATION_COST"));
  params.MAX_SHIFT = (int)mxGetScalar( 
                        mxGetField(params_struct, 0, "MAX_SHIFT"));
  params.NUM_ROTATION = (int)mxGetScalar( 
                        mxGetField(params_struct, 0, "NUM_ROTATION"));
  params.MAX_ROTATION = (float)mxGetScalar( 
                        mxGetField(params_struct, 0, "MAX_ROTATION") );
  params.ROTATION_COST = (float)mxGetScalar( 
                        mxGetField(params_struct, 0, "ROTATION_COST") );    
  params.MAX_STRETCH = (float)mxGetScalar( 
                        mxGetField(params_struct, 0, "MAX_STRETCH") );
  params.STRETCH_COST = (float)mxGetScalar( 
                        mxGetField(params_struct, 0, "STRETCH_COST") );
  params.INVALID_COST = (float)mxGetScalar( 
                        mxGetField(params_struct, 0, "INVALID_COST") );
  params.NUM_STRETCH = 2*params.MAX_STRETCH + 1;
  params.NUM_SHIFT = params.MAX_SHIFT*2+1;

//   mexPrintf("edge-distance-cost: %f\n", params.EDGE_DISTANCE_COST);
//   mexPrintf("translation-cost: %f\n",   params.TRANSLATION_COST);
//   mexPrintf("max-shift: %d\n",          params.MAX_SHIFT);
//   mexPrintf("num-shift: %d\n",          params.NUM_SHIFT);
//   mexPrintf("num-rotation: %d\n",       params.NUM_ROTATION);
//   mexPrintf("max-rotation: %f\n",       params.MAX_ROTATION);
//   mexPrintf("rotation-cost: %f\n",      params.ROTATION_COST);
//   mexPrintf("max-stretch: %f\n",        params.MAX_STRETCH);
//   mexPrintf("num-stretch: %d\n",        params.NUM_STRETCH);
//   mexPrintf("stretch-cost: %f\n",       params.STRETCH_COST);
//   mexPrintf("invalid-cost: %f\n",       params.INVALID_COST);
  
  // --------------------------------------------------------------------//
  sz = (mwSize*)mxGetDimensions( prhs[0] );
  int rows = sz[0];
  int cols = sz[1];
  int chan = sz[2];
  if (chan != 3) { mexPrintf("not a 3-channel image\n"); return; }
  
  plhs[0] = mxCreateDoubleMatrix(ncons, 2, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(ncons, 2, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(ncons, 1, mxREAL);
  
  double* pts_occd_out = (double*)mxGetPr(plhs[0]);
  double* pts_occr_out = (double*)mxGetPr(plhs[1]); 
  double* gmm_weights = (double*)mxGetPr(plhs[2]); 
   
  vector<int> shifts_ = vector<int>(params.NUM_SHIFT);
  for (int i = 0; i < params.NUM_SHIFT; ++i) { shifts_[i] = -params.MAX_SHIFT + i; }
  vector<double> rotations_ = vector<double>(params.NUM_ROTATION);
  for (int i = 0; i < params.NUM_ROTATION; ++i) {
      rotations_[i] = -params.MAX_ROTATION 
                + i*(2*params.MAX_ROTATION/double(params.NUM_ROTATION-1));
  }
  if (params.NUM_ROTATION==1) { rotations_[0] = 0.0; }
   
  // initialize 'invalid' mask
  char* invalid = new char[rows*cols];
  memset(invalid, 0, sizeof(char)*rows*cols );
  
  double occd_pt[2];
  double occr_pt[2];

  for (int i = 0; i < ncons; ++i) { gmm_weights[i] = 1.0; }

  memcpy(pts_occd_out, pts_occd, ncons*2*sizeof(double) );
  memcpy(pts_occr_out, pts_occr, ncons*2*sizeof(double) );

  // how many groups?
  int ngroups = get_num_groups(groups, ncons);

  // get groups assignments?
  vector<vector<int> > cons_group_assignments;
  get_group_assignments(cons_group_assignments, groups, ngroups, ncons);

  //-----------------------------------------------------------------------
  for (int grp = 0; grp < ngroups; ++grp) {       
    vector<int> idx = cons_group_assignments[grp];
    size_t num_pts = idx.size();
    
    double* gmm_fg_mu_  = gmm_fg_mu + grp*K*gmm_dim;
    double* gmm_fg_cov_ = gmm_fg_cov + grp*K*gmm_dim;
    double* gmm_fg_pi_  = gmm_fg_pi + grp*K;

    double* gmm_bg_mu_  = gmm_bg_mu + grp*K*gmm_dim;
    double* gmm_bg_cov_ = gmm_bg_cov + grp*K*gmm_dim;
    double* gmm_bg_pi_  = gmm_bg_pi + grp*K;    
    
    constraint_data data;
    data.gmm_fg.mu  = gmm_fg_mu + grp*K*gmm_dim;
    data.gmm_fg.cov = gmm_fg_cov + grp*K*gmm_dim;
    data.gmm_fg.pi  = gmm_fg_pi + grp*K;
    data.gmm_fg.dim = gmm_dim;
    data.gmm_fg.K = K;
    //
    data.gmm_bg.mu  = gmm_bg_mu + grp*K*gmm_dim;
    data.gmm_bg.cov = gmm_bg_cov + grp*K*gmm_dim;
    data.gmm_bg.pi  = gmm_bg_pi + grp*K;   
    data.gmm_bg.dim = gmm_dim;
    data.gmm_bg.K = K;
    //
    data.occd_pts = vector<double>(2*num_pts);
    data.occr_pts = vector<double>(2*num_pts);
    data.weights = vector<double>(2*num_pts);
    data.separation = vector<double>(num_pts);
    data.num_pts = num_pts;
    
    for (size_t i = 0; i < idx.size(); ++i) { 
        int q = idx[i];
        data.weights[i] = 1.0;
        //
        data.occd_pts[2*i] = MAT2C( pts_occd[ q ] );
        data.occd_pts[2*i+1] = MAT2C( pts_occd[ q+ncons ] );
        data.occr_pts[2*i] = MAT2C( pts_occr[ q ] );
        data.occr_pts[2*i+1] = MAT2C( pts_occr[ q+ncons ] );
        data.separation[i] = max(0.0, pts_max_separation[q]);
    }
    
    if (is_gmm_bad(&data.gmm_bg) || is_gmm_bad(&data.gmm_fg) || (num_pts==0) ) { 
        continue;
    }
    // --------------------------------------------------------------------
    
    // actually perform the optimization over deformation
    optimal_constraint_deformation_A2(I0, edge_image, invalid, &data,
                shifts_, rotations_, &params, rows, cols);
    
    // put the points back...
    for (int i = 0; i < num_pts; ++i) {   
      int q = idx[i];
      pts_occd_out[ q ]         = C2MAT( data.occd_pts[2*i+0] );
      pts_occd_out[ q+ncons ]   = C2MAT( data.occd_pts[2*i+1] );
      pts_occr_out[ q ]         = C2MAT( data.occr_pts[2*i+0] ); 
      pts_occr_out[ q+ncons ]   = C2MAT( data.occr_pts[2*i+1] );
      gmm_weights[ q ] = data.weights[i];
    }
  }
  delete invalid;
   
  //-----------------------------------------------------------------------
}

void optimal_constraint_deformation_A2(double* I0, double* edge_image,
        char* invalid, constraint_data* data,
        vector<int>& shifts_, vector<double>& rotations_, 
        const perturb_params* params, int rows, int cols) {

  // 
  int num_pts = data->num_pts;
  int K = data->gmm_fg.K;
  int dim = data->gmm_fg.dim;
  
  double center[2];
  double occd[2];
  double occr[2];
  double occd_nocenter[2];
  double occr_nocenter[2];
  double occd_nocenter_rot[2];
  double occr_nocenter_rot[2];
  double occd_perturbed[2];
  double occr_perturbed[2];
   
  double I_occd[3];
  double I_occr[3];
  
  double* fg_invsqrtcov = new double[dim*K];
  double* bg_invsqrtcov = new double[dim*K];
  double* fg_logsqrtcov = new double[K];
  double* bg_logsqrtcov = new double[K];
 
  for (int i = 0; i < dim*K; ++i) { 
      fg_invsqrtcov[i] = 1.0/sqrt( data->gmm_fg.cov[i] );
      bg_invsqrtcov[i] = 1.0/sqrt( data->gmm_bg.cov[i] );
  }
  
  for (int k = 0; k < K; ++k) { 
      fg_logsqrtcov[k] = 0;
      bg_logsqrtcov[k] = 0;
      for (int j = 0; j < dim; ++j) { 
          fg_logsqrtcov[k] += 0.5*log( data->gmm_fg.cov[j + k*dim]);
          bg_logsqrtcov[k] += 0.5*log( data->gmm_bg.cov[j + k*dim]);
      }
  }
  
  double p_occr_fg, p_occr_bg, p_occd_fg, p_occd_bg;
  
  xform best_xform;
  best_xform.scale = 1;
  best_xform.theta = 0;
  best_xform.T[0] = 0;
  best_xform.T[1] = 0;

  // ----------------------------------------------------------------------
  // initialize memoization table.
  memorization_table memtable_occd = init_memorization_table(rows, cols, 1);
  memorization_table memtable_occr = init_memorization_table(rows, cols, 1);

  // ----------------------------------------------------------------------
  for (int pt = 0; pt < num_pts; ++pt) { 
      double best_weight = 1.0;
      occd[0] = data->occd_pts[2*pt+0];
      occd[1] = data->occd_pts[2*pt+1];
      occr[0] = data->occr_pts[2*pt+0];
      occr[1] = data->occr_pts[2*pt+1];
      // compute center
      center[0] = (occd[0]+occr[0])/2.0;
      center[1] = (occd[1]+occr[1])/2.0;
      // remove center:
      occd_nocenter[0] = occd[0]-center[0];
      occd_nocenter[1] = occd[1]-center[1];
      occr_nocenter[0] = occr[0]-center[0];
      occr_nocenter[1] = occr[1]-center[1];
      double dx = (occd[0]-occr[0]);
      double dy = (occd[1]-occr[1]);
      double separation = sqrt(dx*dx + dy*dy);

      double Emaxval = -99999999;
      
      // ------------------------------------------------------------------    
      // this block ensures that it is possible to evaluate stretch
      // which is smaller than max-allowed-separation.
      // ------------------------------------------------------------------
      int NUM_STRETCH = params->NUM_STRETCH;
      int MAX_STRETCH = params->MAX_STRETCH;
      adjust_num_stretch(&NUM_STRETCH, &MAX_STRETCH, 
                            separation, data->separation[pt]);                      
      // ------------------------------------------------------------------
      for (int tt = 0; tt < NUM_STRETCH; ++tt) { 
          int stretch = tt - MAX_STRETCH;
          double scale = max(0.0, (separation+stretch)/(separation) );
          // if distance between the points becomes too large, skip the step..
          if (separation < 1e-5) { continue; }
          if (scale == 0) { continue; }
          if (mxIsNaN(scale)) { continue; }
          if (mxIsInf(abs(scale))) { continue; }

          // it is also pointless to check this case.
          if (scale*separation > data->separation[pt]) { continue; }

          for (int k = 0; k < rotations_.size(); ++k) {
              double theta = rotations_[k];
              rotate_and_scale( occd_nocenter_rot, occd_nocenter, center, theta, scale );
              rotate_and_scale( occr_nocenter_rot, occr_nocenter, center, theta, scale );
              //
              int x1_,y1_,x2_,y2_;
              for (int ii = 0; ii < shifts_.size(); ++ii) { 
                x1_ = occd_nocenter_rot[0] + shifts_[ii];
                x2_ = occr_nocenter_rot[0] + shifts_[ii];
                x1_ = min( cols-1, max( x1_, 0) );
                x2_ = min( cols-1, max( x2_, 0) );
 
                for (int jj = 0; jj < shifts_.size(); ++jj) { 
                    y1_ = occd_nocenter_rot[1] + shifts_[jj];
                    y2_ = occr_nocenter_rot[1] + shifts_[jj];
                    y1_ = min( rows-1, max(y1_, 0) );
                    y2_ = min( rows-1, max(y2_, 0) );
                    
                    double p_occd, p_occr, prob_occr, prob_occd;
                    short memtable_idx;

                    int lin_idx1 = linear_index(y1_, x1_, rows); // occluder
                    int lin_idx2 = linear_index(y2_, x2_, rows); // occluder
                    double d1 = edge_image[ lin_idx1 ];
                    double d2 = edge_image[ lin_idx2 ];
                    
                    double invalid_cost = 0;
                    if ((invalid[ lin_idx1 ]==1) || (invalid[ lin_idx2 ]==1)) {
                        invalid_cost = params->INVALID_COST;
                    }

                    I_occd[0] = I0[ linear_index(y1_, x1_, 0, rows, cols) ];
                    I_occr[0] = I0[ linear_index(y2_, x2_, 0, rows, cols) ];         
                    I_occd[1] = I0[ linear_index(y1_, x1_, 1, rows, cols) ];
                    I_occr[1] = I0[ linear_index(y2_, x2_, 1, rows, cols) ];         
                    I_occd[2] = I0[ linear_index(y1_, x1_, 2, rows, cols) ];
                    I_occr[2] = I0[ linear_index(y2_, x2_, 2, rows, cols) ];         

                    memtable_idx = is_in_memorization_table(memtable_occd, y1_, x1_, 0); // is pt=0 ?
                    if (memtable_idx>0) {
                        p_occd = memtable_occd.data_log[memtable_idx];
                        prob_occd = memtable_occd.data_prob[memtable_idx];
                    } else {
                        // occluded - foreground
                        p_occd_fg = utils_eval_gmm_fast( I_occd, 
                                         data->gmm_fg.mu, data->gmm_fg.cov,
                                         fg_invsqrtcov, fg_logsqrtcov, 
                                         data->gmm_fg.pi, dim, K);
                        // occluded - background
                        p_occd_bg = utils_eval_gmm_fast( I_occd, 
                                         data->gmm_bg.mu, data->gmm_bg.cov, 
                                         bg_invsqrtcov, bg_logsqrtcov, 
                                         data->gmm_bg.pi, dim, K);

                        p_occd = convert_to_log_ratio(p_occd_bg, p_occd_fg);
                        prob_occd = p_occd_bg / (p_occd_fg + p_occd_bg);

                        put_in_memorization_table(memtable_occd, p_occd, prob_occd, y1_, x1_, 0); // is pt=0 ?
                    }                    
                    
                    memtable_idx = is_in_memorization_table(memtable_occr, y2_, x2_, 0); // is pt=0 ?               
                    
                    if (memtable_idx>0) {
                        p_occr = memtable_occr.data_log[memtable_idx];
                        prob_occr =  memtable_occr.data_prob[memtable_idx];                      
                    } else {
                        // occluder - foreground
                        p_occr_fg = utils_eval_gmm_fast( I_occr, 
                                     data->gmm_fg.mu, data->gmm_fg.cov, 
                                     fg_invsqrtcov, fg_logsqrtcov, 
                                     data->gmm_fg.pi, dim, K);
                        // occluder - background
                        p_occr_bg = utils_eval_gmm_fast( I_occr, 
                                     data->gmm_bg.mu, data->gmm_bg.cov, 
                                     bg_invsqrtcov, bg_logsqrtcov, 
                                     data->gmm_bg.pi, dim, K);

                        p_occr = convert_to_log_ratio(p_occr_fg, p_occr_bg);
                        prob_occr = p_occr_fg / (p_occr_fg+p_occr_bg);

                        put_in_memorization_table(memtable_occr, p_occr, prob_occr, y2_, x2_, 0); // is pt=0 ?
                    }

                    double dist = shifts_[jj]*shifts_[jj] + shifts_[ii]*shifts_[ii];
                    double deformation_cost = 
                                 params->TRANSLATION_COST*dist + 
                                 params->ROTATION_COST*abs(theta) + 
                                 params->STRETCH_COST*abs(stretch);

                    double energy = ( (p_occr+p_occd) - deformation_cost 
                                 - invalid_cost - params->EDGE_DISTANCE_COST*(d1+d2) );

                    // ---------------------------------------------------
                    // store best transformation..
                    if (energy > Emaxval) { 
                        best_xform.theta = theta;
                        best_xform.scale = scale;
                        best_xform.T[0] = shifts_[ii];
                        best_xform.T[1] = shifts_[jj];
                        best_weight = compute_gmm_weight(prob_occr, prob_occd);
                        Emaxval = energy;
                    }
                    // ---------------------------------------------------      
                  }
              }
          }
    }
    //.....................................................................
    // apply the best transform to points:
    double temp[2];
    memset(temp, 0, sizeof(double)*2);
    rotate_and_scale( temp, occd_nocenter, center, 
                                     best_xform.theta, best_xform.scale );
    data->occd_pts[2*pt+0] = floor( temp[0]+best_xform.T[0] );
    data->occd_pts[2*pt+1] = floor( temp[1]+best_xform.T[1] );
    rotate_and_scale( temp, occr_nocenter, center, 
                                    best_xform.theta, best_xform.scale );
    data->occr_pts[2*pt+0] = floor( temp[0]+best_xform.T[0] );
    data->occr_pts[2*pt+1] = floor( temp[1]+best_xform.T[1] );
    
    force_inside(&data->occd_pts[2*pt], rows, cols);
    force_inside(&data->occr_pts[2*pt], rows, cols);    
    data->weights[pt] = best_weight;
    
    // update invalid mask
    invalid[ linear_index( data->occd_pts[2*pt+1], data->occd_pts[2*pt+0], rows) ] = 1;
    invalid[ linear_index( data->occr_pts[2*pt+1], data->occr_pts[2*pt+0], rows) ] = 1;
  }
  
  //
  delete fg_invsqrtcov;
  delete bg_invsqrtcov;
  delete fg_logsqrtcov;
  delete bg_logsqrtcov; 
}

void force_inside(double* pt, int rows, int cols) { 
    pt[0] = min(max(0.0 , pt[0]), (double)(cols-1) );
    pt[1] = min(max(0.0 , pt[1]), (double)(rows-1) );
}


void rotate_and_scale( double* out, double* in, double* center, double theta, double scale ) {
    out[0] = (cos(+theta)*in[0] + sin(+theta)*in[1])*scale + center[0];
    out[1] = (sin(-theta)*in[0] + cos(+theta)*in[1])*scale + center[1];    
}

memorization_table init_memorization_table(int rows, int cols, int num_pts) {
    memorization_table m;
    m.n = num_pts;
    m.rows = rows;
    m.cols = cols;

    // stores indices into arrays below:
    m.idx = vector<short>(num_pts*rows*cols);
    memset(&m.idx[0], 0, sizeof(short)*num_pts*rows*cols);
    m.num = 0;
    m.data_log= vector<double>();
    m.data_prob = vector<double>();
    m.data_log.reserve(256);
    m.data_prob.reserve(256);
    // dummy, since counting starts at 1.
    m.data_log.push_back( - 1 );
    m.data_prob.push_back( - 1 );
    
    return m;
}

short inline is_in_memorization_table(memorization_table& m, int r, int c, int pt) {
    return m.idx[ linear_index(r, c, pt, m.rows, m.cols) ];
}

void inline put_in_memorization_table(memorization_table& m, double log, double prob, int r, int c, int pt) { 
    m.num++;
    m.idx[ linear_index(r,c,pt,m.rows,m.cols) ] = m.num;
    m.data_log.push_back( log );
    m.data_prob.push_back( prob );
}

int get_num_groups(const double* groups, int n) {
  int ngroups = 0;
  for (int i = 0; i < n; ++i) { ngroups = max((int)groups[i],ngroups); }
  return ngroups;
}

void get_group_assignments(vector<vector<int> >& assignments, double* groups, int ngroups, int n) { 
    int k = 1;
    assignments = vector<vector<int> >(ngroups);
    if ((n>0) && (ngroups>0)) { k = n/ngroups + 1; }

    for (int i = 0; i < ngroups; ++i) { 
        assignments[i].reserve(k);
    }    
    // --------------------------------
    for (int i = 0; i < n; ++i) { 
        int grp = MAT2C( groups[i] );
        assignments[grp].push_back(i);
    }
}

int is_gmm_bad(gmm* g) { 
    int bad1 = (g->mu == 0) || (mxIsNaN(g->mu[0]) || mxIsNaN(g->mu[0]));
    int bad2 = (g->cov == 0) || (mxIsNaN(g->cov[0]) || mxIsNaN(g->cov[0]));
    int bad3 = (g->pi == 0) || (mxIsNaN(g->pi[0]) || mxIsNaN(g->pi[0]));
    return (bad1 || bad2 || bad3);
}

void adjust_num_stretch(int* NUM_STRETCH, int* MAX_STRETCH, 
                            double separation, double maxallowedseparation) { 
    double min_separation = get_min_nonzero_separation(separation, *NUM_STRETCH, *MAX_STRETCH);
    while (min_separation > maxallowedseparation) {
      *MAX_STRETCH += 1;
      *NUM_STRETCH = 2*(*MAX_STRETCH) + 1;
      min_separation = get_min_nonzero_separation(separation, *NUM_STRETCH, *MAX_STRETCH);
    }
}

double get_min_nonzero_separation(double sep, int NUM_STRETCH, int MAX_STRETCH) {
    double min_nonzero_separation = 0;
    for (int tt = 0; tt < NUM_STRETCH; ++tt) { 
       min_nonzero_separation = (sep-MAX_STRETCH+tt);
       if (min_nonzero_separation>0) { break; }
    }
    return min_nonzero_separation;
}
