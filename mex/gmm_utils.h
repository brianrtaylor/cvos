#ifndef GMM_UTILS_H_
#define GMM_UTILS_H_

#include <cstdlib>
#include <cmath>
#include <algorithm> // max/min

#define LOG2PI 1.837877066409345
#define PI 3.141592653589793
using namespace std;

double compute_gmm_weight(double p_occr_fg, double p_occr_bg, double p_occd_fg, double p_occd_bg);

double inline compute_gmm_weight(double P_occr, double P_occd) { 
    double p = (P_occr+P_occd)/2.0;
    return 2.0*max(0.0, p-0.5);    
}

double convert_to_log_ratio(double p1, double p2);

double normpdf(double* I, double* mu, double* cov, int dim);

double utils_eval_gmm(double* I, double* mu, double* cov, double* pi, int dim, int K);

double compute_gmm_weight(double p_occr_fg, double p_occr_bg, double p_occd_fg, double p_occd_bg);

double utils_eval_gmm_fast(double* I, double* mu, double* cov, double* pi, int dim, int K);

double utils_eval_gmm_fast(const double* I, const double* mu, 
        const double* cov, const double* invsqrtcov, 
        const double* logsqrtcov, const double* pi, int dim, int K);

#endif // GMM_UTILS_H_