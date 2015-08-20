// @file tvsegment.cpp
#include "tvsegment.h"

tvsegment::tvsegment() { 
}

tvsegment::~tvsegment() {
}

tvsegment::tvsegment(const vector<float>& img, tvsegment_params* p) { 
    rows = p->rows;
    cols = p->cols;
    chan = p->chan;
    K = p->num_clusters;
    
    // sigma_c - set this
    // sigma_y - set this
    float sigma_u = 0.95f/sqrt(8);
    float sigma_y = 0.95f/sigma_u/8.0f;

    // compute tv regularization weights and weigh it by lambda 
    weights = compute_weights(img, p->beta);
    for (size_t i = 0; i < weights.size(); ++i) { weights[i] *= p->lambda; }
    // vis_weights(weights);
    
    // initialize optimization variables
    u = vector<float>(rows*cols*K, 0); // nodes
    vector<float> y1 = vector<float>(rows*cols*K*2, 0); // edges
    vector<float> x = vector<float>(rows*cols*K, 0); // dual variable
    
    vector<float> dx = vector<float>(rows*cols*K*2, 0); // nabla u
    vector<float> dxt = vector<float>(rows*cols*K, 0); // div 
    vector<float> wnorm;

    if (p->isotropic) { 
        wnorm = vector<float>(rows*cols*K);
        int m = rows*cols*K;
        for (int i = 0; i < rows*cols*K; ++i) { 
            wnorm[i] = sqrt(weights[i]*weights[i] + weights[i+m]*weights[i+m]);
        }
    }
    for (int t = 0; t < p->max_iters; ++t) { 
//         if ((t%10)==0) printf("%d/%d\n", t, p->max_iters );
        // apply pixelwise gradient operator:
        apply_pixelwise_gradient_op(dx, x, rows, cols, K);
        if (p->isotropic == 0) { // prox operator for anisotropic TV penalty
            for (size_t i = 0; i < y1.size(); ++i) { 
                y1[i] = y1[i] + dx[i]*sigma_y;
                y1[i] = -BCV_SIGN(y1[i])*max( -weights[i], -abs(y1[i]) );
            }
        } else { // isotropic TV penalty:
            int m = rows*cols*K;
            for (size_t i = 0; i < y1.size()/2; ++i) { 
                y1[i] = y1[i] + dx[i]*sigma_y;
                y1[i+m] = y1[i+m] + dx[i+m]*sigma_y;
                float q1 = y1[i];    
                float q2 = y1[i+m];
                float qn = sqrt(q1*q1+q2*q2);
                float f = min( 1.0f, wnorm[i]/qn ); 
                y1[i] = q1*f;
                y1[i+m] = q2*f;
            }
        }
        // 
        apply_pixelwise_gradient_op_transpose(dxt, y1, rows, cols, K);
        for (size_t i = 0; i < u.size(); ++i) {
            x[i] = u[i];
            u[i] = u[i] - sigma_u*(dxt[i] + p->unary[i]);
        }
        project_onto_prob_simplex(u);
        for (size_t i = 0; i < u.size(); ++i) { 
            x[i] = 2*u[i] - x[i];   
        }
        //vis_assignments(u);
    }
}

vector<int> tvsegment::get_assignments() {
    vector<int> seg = vector<int>(rows*cols, -1);
    int idx;
    for (int i = 0; i < seg.size(); ++i) { 
        //
        float maxval = -1.0f;
        
        for (int k = 0; k < K; ++k) { 
            #if defined(HAVE_MATLAB)
            idx = i+rows*cols*k;
            #else
            idx = K*i+k;
            #endif       
            if (u[idx] > maxval) {
                maxval = u[idx];
                seg[i] = k;
            }
        }
    }
    return seg;
}

vector<float> tvsegment::get_weights() { 
    return weights;
}
    
// computes exp(-|I(x)-I(y)|^2/beta), with beta = E( |I(x)-I(y)|^2 )
vector<float> tvsegment::compute_weights(const vector<float>& img, float beta) { 
    vector<float> weights = vector<float>(rows*cols*K*2); 

    float inv_den_sqrt = 1.0f/sqrt(float(rows*cols*K*2));
    
    vector<float> grad = vector<float>(rows*cols*chan*2);
    apply_pixelwise_gradient_op(grad, img, rows, cols, chan);
    float Z = 0.0f;
    for (size_t i = 0; i < grad.size(); ++i) { 
        Z += grad[i]*grad[i]*inv_den_sqrt;
    }
    Z *= inv_den_sqrt;
    Z *= beta;
    // Z = beta/n \sum_i grad(i)^2
    
    for (int d = 0; d < 2; ++d) { 
        int diff_offset_in = d*rows*cols*chan;
        int diff_offset_out = d*rows*cols*K;
        
        for (int i = 0; i < rows*cols; ++i) { 
            float dist = 0;
            for (int k = 0; k < chan; ++k) {
                #if defined(HAVE_MATLAB)
                int idx = i+rows*cols*k + diff_offset_in; // index into 'grad'
                #else
                int idx = chan*i+k + diff_offset_in;
                #endif
                dist += grad[idx]*grad[idx];
            }
            for (int k = 0; k < K; ++k) {
                #if defined(HAVE_MATLAB)
                int idx = i + rows*cols*k + diff_offset_out;
                #else
                int idx = K*i + k + diff_offset_out;
                #endif
                weights[idx] = exp(-0.5f*dist/Z);
            }
        }
    }
    return weights;
}

void tvsegment::project_onto_prob_simplex(vector<float>& x) {
    #ifdef HAVE_MATLAB
    project_onto_prob_simplex_matlab(x);
    return;
    #else
    int n = x.size()/K;
    float* p;
    for (int i = 0; i < n; ++i) { // for each pixel
        // project onto simplex:
        p = (float*)&x[K*i];
        float mu = -1.0f;
        for (int k1 = 0; k1 < K; ++k1) { 
            float a = 0;
            int b = 0;
            for (int k2 = 0; k2 < K; ++k2) {
                if (p[k2] > mu) {
                    a += p[k2];
                    b += 1;
                }
            }
            mu = (a-1)/b;
        }
        for (size_t k = 0; k < K; ++k) { 
            p[k] = max(p[k] - mu, 0.0f);
        }
    }
    #endif
}

void tvsegment::project_onto_prob_simplex_matlab(vector<float>& x) {
    int n = x.size()/K;
    vector<float> p = vector<float>(K);
    for (int i = 0; i < n; ++i) { // for each pixel
        // project onto simplex:
        for (int k = 0; k < K; ++k) { p[k] = x[i+rows*cols*k]; }
        float mu = -1.0f;
        for (int k1 = 0; k1 < K; ++k1) { 
            float a = 0;
            int b = 0;
            for (int k2 = 0; k2 < K; ++k2) {
                if (p[k2] > mu) {
                    a += p[k2];
                    b += 1;
                }
            }
            mu = (a-1)/b;
        }
        for (size_t k = 0; k < K; ++k) { 
            p[k] = max(p[k] - mu, 0.0f);
            x[i+rows*cols*k] = p[k];
        }
    }
}

#ifdef HAVE_OPENCV
void tvsegment::vis_weights(const vector<float>& weights) { 
    vector<float> img = vector<float>(rows*cols);
    for (int i = 0; i < img.size(); ++i) { 
        img[i] = (weights[K*i] + weights[rows*cols*K + K*i])/2.0f;
    }
    bcv_imshow<float>("weights", img, rows, cols, 1, CV_32F);
    cv::waitKey(0);
}

void tvsegment::vis_assignments(const vector<float>& u) { 
    vector<float> img = vector<float>(rows*cols*K);
    int i = 0;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) { 
            for (int k = 0; k < K; ++k) { 
                img[ linear_index(r,c+cols*k,cols*K) ] = u[i];
                i++;
            }
        }
    }
    bcv_imshow<float>("labels", img, rows, cols*K, 1, CV_32F);
    cv::waitKey(0);
}
#endif
