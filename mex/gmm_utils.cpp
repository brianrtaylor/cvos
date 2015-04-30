#include "gmm_utils.h"

double compute_gmm_weight(double p_occr_fg, double p_occr_bg, double p_occd_fg, double p_occd_bg) { 
    double P_occr = p_occr_fg / (p_occr_fg + p_occr_bg);
    double P_occd = p_occd_bg / (p_occd_fg + p_occd_bg);
    double p = (P_occr+P_occd)/2.0;    
    return 2.0*max(0.0, p-0.5);
}

double convert_to_log_ratio(double p1, double p2) {
    p1 = max(min(300.0, log(p1)), -300.0);
    p2 = max(min(300.0, log(p2)), -300.0);
    return (p1-p2);
}

double utils_eval_gmm_fast(const double* I, const double* mu, 
        const double* cov, const double* invsqrtcov, 
        const double* logsqrtcov, const double* pi, int dim, int K) {
    double p = 0;
    double f;
    const double *MU, *COV, *INVSQRTCOV;
    for (int k = 0; k < K; ++k) {
        MU = mu + dim*k;
        COV = cov + dim*k;     
        INVSQRTCOV = invsqrtcov + dim*k;
        
        double logp = -dim*LOG2PI - logsqrtcov[k];

        for (int i = 0; i < dim; ++i) { 
            f = (I[i]-MU[i])*INVSQRTCOV[i];
            logp -= 0.5*f*f;
        }
        p += pi[k]*exp(logp);
    }
    return p;
}


double utils_eval_gmm_fast(double* I, double* mu, double* cov, double* pi, int dim, int K) {
  double p = 1e-10;
  double *MU, *COV;
  for (int k = 0; k < K; ++k) {
    MU = mu + dim*k;
    COV = cov + dim*k;     
    double kp = 1.0;
    for (int i = 0; i < dim; ++i) { 
      kp *= (1.0 / sqrt(COV[i] * 2 * PI));
      kp *= exp((-0.5 * (I[i] - MU[i])*(I[i] - MU[i])) / COV[i]);
    }
    p += pi[k]*kp;
  }
  return p;
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
    double logp = -dim*LOG2PI; // factor of 2 is unnecessary
    double sqrtcov, f;
    for (int i = 0; i < dim; ++i) { 
        sqrtcov = sqrt(cov[i]);
        f = (I[i]-mu[i])/sqrtcov;
        logp += (-0.5*f*f - log(sqrtcov));
    }
    return exp(logp);
}