#ifndef BCV_TVSEGMENT_H_
#define BCV_TVSEGMENT_H_

#include <cstdlib>
#include <cmath>
#include <limits>
#include "bcv_diff_ops.h"
#ifdef HAVE_OPENCV
#include "bcv_opencv.cpp"
using namespace cv;
#endif

#ifndef BCV_SIGN
#define BCV_SIGN(x) (((x) > 0) - ((x) < 0)) // -1, 0, or +1
#endif

using namespace std;

//! tv segmentation parameters.
struct tvsegment_params {
    float lambda; //! weight of TV penalty
    float beta; //! sigma of weights. small value creates more extreme weights
    vector<float> unary; // unary term
    vector<float> clusters; // values of clusters (only for visualization)
    int rows;
    int cols;
    int chan; // this is also the dimension of each cluster
    int num_clusters;
    int max_iters; //! maximum number of iterations
    int isotropic; //! perform isotropic TV regularization or not
};

class tvsegment {
    public:
        //! main function that performs optimization
        tvsegment(const vector<float>& img, tvsegment_params* p);
        ~tvsegment();
        tvsegment();
        vector<int> get_assignments();
        vector<float> get_weights();
    private:
        vector<float> u;
        vector<float> weights;
        int rows;
        int cols;
        int chan;
        int K;

        //! computes weights for tv regularization
        vector<float> compute_weights(const vector<float>& img, float beta);
 
        void project_onto_prob_simplex(vector<float>& x);
        #ifdef HAVE_MATLAB
        void project_onto_prob_simplex_matlab(vector<float>& x);
        #endif

        // visualization (mostly for debugging)
        #ifdef HAVE_OPENCV
        void vis_weights(const vector<float>& weights);
        void vis_assignments(const vector<float>& u);
        #endif
};

#endif // BCV_TVSEGMENT_H_
