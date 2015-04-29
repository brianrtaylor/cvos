#include "solver_primaldual.h"

solver_primaldual::solver_primaldual() { 
    params = NULL;
    sol_x.clear();
    sol_c.clear();    
}

solver_primaldual::~solver_primaldual() { }

solver_primaldual::solver_primaldual(primaldual_params* p) {
    // stub
    params = p;
}

vector<float> solver_primaldual::solve() {
    if (!check_params()) {
        return vector<float>();
    }
    // initialize variables:
    sol_y1 = vector<float>(params->prob.nedges, 0.0f);
    sol_y2 = vector<float>(params->prob.nocc_constraints, 0.0f);
    if (params->opt.use_temporal_penalty && 
            (params->opt.use_ell_infinity || params->opt.use_both_penalties) ) {
        sol_y3 = vector<float>(params->prob.nnodes, 0.0f);
    }
    
    // only initialize if dimension does not match current dimension
    if (sol_x.size() != (size_t)params->prob.nnodes) {
        sol_x = vector<float>(params->prob.nnodes, 0.0f);
    }
    if (sol_c.size() != (size_t)params->prob.nnodes) { 
        sol_c = vector<float>(params->prob.nnodes, 0.0f);
    }
    if (params->opt.use_initial_layers) { 
        memcpy( &sol_c[0], &params->opt.init_layers[0], sizeof(float)*params->prob.nnodes ); 
        memcpy( &sol_x[0], &params->opt.init_layers[0], sizeof(float)*params->prob.nnodes ); 
    }
    
    // initialize simplex-projection variables (even if L-infinity is not used).
    simplex_proj_lb = 0.0f;
    simplex_proj_ub = 1.0f;
    simplex_proj_sol = 0.0f;
    simplex_scaling = 1.0f;
    
    // temporary (intermediate) variables
    vector<float> DX = vector<float>(params->prob.nedges, 0.0f);
    vector<float> DOCCX = vector<float>(params->prob.nocc_constraints, 0.0f);
    vector<float> DTY = vector<float>(params->prob.nnodes, 0.0f);
    vector<float> DOCCTY = vector<float>(params->prob.nnodes, 0.0f);
    
    float fx=numeric_limits<float>::max();
    float fx_prev, dfx=1, elapsed=0, elapsed_y1=0, elapsed_y2=0, elapsed_c=0;
    int DFX_SIZE = 3;
    float dx_ = 0.0f, dx = 1.0f;
    float fcalc_elapsed = 0;
    vector<float> dfx_vec = vector<float>(DFX_SIZE, numeric_limits<float>::max() );
    vector<float> dx_vec = vector<float>(DFX_SIZE, numeric_limits<float>::max() );
    double t1=0, t2=0;
    float sigma_y = params->opt.sigma_y;
    float sigma_c = params->opt.sigma_c;

    
    for (int iter = 0; iter < params->opt.max_iterations; ++iter) { 
        if (iter % params->opt.verbosity == 0)
            t1 = now_ms();
        // apply difference operator to x. exploiting lattice structure helps a lot
        if (params->prob.solve_pixelwise) {
            apply_pixelwise_gradient_op(DX, sol_x, params->prob.rows, params->prob.cols);
        } else {
            apply_sparse_op<int>(DX, params->prob.D, sol_x );
        }

        apply_sparse_op<int>(DOCCX, params->prob.Docc, sol_x);
        // --------------------------------------------------------------------
        //                        1. solve for y
        // --------------------------------------------------------------------
        // update 1:
        if (iter % params->opt.verbosity == 0) { t2 = now_ms(); }        
        solve_y1(DX, sigma_y);
        if (iter % params->opt.verbosity == 0) { elapsed_y1 = (now_ms()-t2); }
        // update 2:
        if (iter % params->opt.verbosity == 0) { t2 = now_ms(); }
        for (size_t i = 0; i < sol_y2.size(); ++i) {
            float y2_hat = sol_y2[i] + sigma_y * DOCCX[i];
            sol_y2[i] = min( max( y2_hat - sigma_y, -params->prob.occweights[i] ), 0.0f);
        }
        if (iter % params->opt.verbosity == 0) { elapsed_y2 = (now_ms()-t2); }

        // update y3 if necessary:
        if (params->opt.use_temporal_penalty && 
                (params->opt.use_ell_infinity || params->opt.use_both_penalties)) {
            solve_y3(sigma_y);
        }
        
        // --------------------------------------------------------------------
        //                        2. solve for c
        // --------------------------------------------------------------------
        // apply transposed difference operators to y:
        if (params->prob.solve_pixelwise) {
            apply_pixelwise_gradient_op_transpose(
                    DTY, sol_y1, params->prob.rows, params->prob.cols);
        } else {
            apply_sparse_op<int>(DTY, params->prob.D, sol_y1, 'T');     
        }
        // --------------------------------------------------------------------
        apply_sparse_op<int>(DOCCTY, params->prob.Docc, sol_y2, 'T');


        if (iter % params->opt.verbosity == 0) { t2 = now_ms(); }

        if (params->opt.use_both_penalties) {
            // note: this works with and without temporal-penalty (i.e. FG prior) 
            dx_ = solve_c_ell1_ell_infinity(DTY, DOCCTY, sigma_c);     
        } else { 
            // use EITHER L1 OR L-INF PENALTY
            if (params->opt.use_ell_infinity) { 
                // use Linfinity regularization on the layer values. this is SLOW.
                dx_ = solve_c_ell_infinity(DTY, DOCCTY, sigma_c);
            } else {
                // use L1 regularization on the layer values. this is fast.
                if (params->opt.use_temporal_penalty) { 
                    dx_ = solve_c_temporal(DTY, DOCCTY, sigma_c, params->prob.tau1);
                } else {
                    dx_ = solve_c(DTY, DOCCTY, sigma_c, params->prob.tau1);
                }
            }
        }
        if (iter % params->opt.verbosity == 0) {
            elapsed_c = (now_ms()-t2);
            dx_vec[iter % DFX_SIZE] = sqrt(dx_)/float(sol_c.size());
            dx = accumulate(dx_vec.begin(), dx_vec.end(), 0.0f)/float(DFX_SIZE); 
            elapsed = (now_ms()-t1);
        }
        // --------------------------------------------------------------------
        //                 performance diagnostics and stopping
        // --------------------------------------------------------------------


        // function value is calculated only if value change is used as a 
        // stopping criterion. this is done for computational reasons; 
        // function evaluations are fairly expensive.
        fx_prev = fx;
        if (params->opt.fx_tolerance > 0) {
            t2 = now_ms();
            fx = eval_opt_func_value();
            fcalc_elapsed = (now_ms()-t2);
            dfx_vec[ iter % DFX_SIZE ] = abs(fx-fx_prev)/fx;
            dfx = accumulate(dfx_vec.begin(), dfx_vec.end(), 0.0f)/float(DFX_SIZE); 
        } else {
            fx = -1;
            dfx = 1.0;
        }
        if ((params->opt.verbosity>0) && (iter % params->opt.verbosity == 0)) { 
            printf("%4d: fx=%.4e | dfx=%.4e | dx=%.4e | %.3f ms | %.3f ms ||", 
                    iter, fx, dfx, dx, elapsed, fcalc_elapsed);
            printf(" %.3f | %.3f | %.3f |", elapsed_y1, elapsed_y2, elapsed_c);
            if (params->opt.use_ell_infinity || params->opt.use_both_penalties) {
                printf(" %2d ", simplex_num_iters);
            }
            printf("\n");
        }
        if ((iter > 10) && (dfx < params->opt.fx_tolerance) && (dfx==dfx)) {
            if (params->opt.verbosity > 0) {
                printf("Change in f(x) too small: %g < %g\n", 
                        dfx, params->opt.fx_tolerance);
            }
            break;
        }
        if ((iter > 10) && (dx < params->opt.dx_tolerance)) {
            if (params->opt.verbosity > 0) {
                printf("Change in x too small: %g < %g\n", 
                        dx, params->opt.dx_tolerance);
            }
            break;
        }
        if (iter >= params->opt.max_iterations) {
            if (params->opt.verbosity > 0) {
                printf("Reached maximum number of iterations.\n");
            }
            break;
        }

        #ifdef HAVE_MATLAB
        if ((params->opt.fx_truth !=0) && (params->opt.fx_truth_eps !=0)) {
            float gap = abs(fx - params->opt.fx_truth)/params->opt.fx_truth;
            if (gap < params->opt.fx_truth_eps) {
                break;
            }
        }
        #endif
    }
    return sol_c;
}

float solver_primaldual::eval_opt_func_value() {
    vector<float> Doccx;
    vector<float> Dx;
    if (params->prob.solve_pixelwise) {
        apply_pixelwise_gradient_op(Dx, sol_c, params->prob.rows, params->prob.cols);
    } else {
        apply_sparse_op<int>(Dx, params->prob.D, sol_c);
    }
    apply_sparse_op<int>(Doccx, params->prob.Docc, sol_c);
    //
    float fx = 0.0f;
    if (params->opt.use_both_penalties) { 
        // \tau 1^T c
        for (size_t i = 0; i < sol_c.size(); ++i) { fx += params->prob.tau1[i] * sol_c[i]; }
        // \tau max(c)
        fx += params->prob.tau2 * (*max_element( sol_c.begin(), sol_c.end() ));
    } else {
        if (params->opt.use_ell_infinity) {
            // \tau max(c)
            fx += params->prob.tau2 * (*max_element( sol_c.begin(), sol_c.end() ));
        } else {
            // \tau 1^T c
            for (size_t i = 0; i < sol_c.size(); ++i) { fx += params->prob.tau1[i] * sol_c[i]; }
        }
    }
    // \lambda^T max(0, 1 - Docc x)
    for (size_t i = 0; i < Doccx.size(); ++i) { 
        fx += params->prob.occweights[i] * max( 0.0f, 1.0f - Doccx[i] );
    }
    // |WDx|_1
    for (size_t i = 0; i < Dx.size(); ++i) {
        fx += params->prob.weights[i] * abs(Dx[i]);
    }
    // kappa^T max(0, 1-c)
    if (params->opt.use_temporal_penalty) { 
        for (size_t i = 0; i < params->prob.kappa.size(); ++i) { 
            fx += params->prob.kappa[i] * max( 0.0f, 1.0f - sol_c[i] );
        }
    }
    return fx;
}

void solver_primaldual::print_cost() {
    //
    float cumfx = 0.0f;
    float fx = 0.0f;

    vector<float> Dx;
    vector<float> Doccx;

    if (params->prob.solve_pixelwise) {
        apply_pixelwise_gradient_op(Dx, sol_c, params->prob.rows, params->prob.cols);
    } else {
        apply_sparse_op<int>(Dx, params->prob.D, sol_c );
    }
    apply_sparse_op<int>(Doccx, params->prob.Docc, sol_c);
        
    if (params->opt.use_both_penalties) { 
        // \tau 1^T c
        fx = 0.0; for (size_t i = 0; i < sol_c.size(); ++i) { fx += params->prob.tau1[i] * sol_c[i]; }
        printf("tau1^T c \t\t\t= %f\n", fx);
        cumfx += fx;
        // \tau 1^T c
        fx = params->prob.tau2 * (*max_element( sol_c.begin(), sol_c.end() ) );
        printf("tau2 max(c) \t\t\t= %f\n", fx);
        cumfx += fx;
    } else {
        if (params->opt.use_ell_infinity) { 
            // \tau 1^T c
            fx = params->prob.tau2 * (*max_element( sol_c.begin(), sol_c.end() ) );
            cumfx += fx;
            printf("tau2 max(c) \t\t\t= %f\n", fx);            
        } else {
            // \tau 1^T c
            fx = 0; for (size_t i = 0; i < sol_c.size(); ++i) { fx += params->prob.tau1[i] * sol_c[i]; }
            cumfx += fx;
            printf("tau1^T c \t\t\t= %f\n", fx);
        }
    }
    // \lambda^T max(0, 1 - Docc x)
    fx = 0;
    for (size_t i = 0; i < Doccx.size(); ++i) { 
        fx += params->prob.occweights[i] * max( 0.0f, 1.0f - Doccx[i] );
    }
    printf("lambda^T max(0, 1-Docc c) \t= %f\n", fx);
    cumfx += fx;
    
    // |WDx|_1
    fx = 0;
    for (size_t i = 0; i < Dx.size(); ++i) {
        fx += params->prob.weights[i] * abs(Dx[i]);
    }
    printf("|| W D c ||_1 \t\t\t= %f\n", fx);
    cumfx += fx;

    // kappa^T max(0,1-c)
    if (params->opt.use_temporal_penalty) { 
        fx = 0;
        for (size_t i = 0; i < sol_c.size(); ++i) {
            fx += params->prob.kappa[i] * max(0.0f, 1-sol_c[i]);
        }
        printf("kappa^T max(0,1-c) \t\t= %f\n", fx);
        cumfx += fx;
    }
    printf("total \t\t\t\t= %f\n", cumfx);
}


void solver_primaldual::project_onto_prob_simplex(vector<float>& x, float eta, float tol) {
    float base_scale = 1.025;
    if (simplex_num_iters > 10) {
        simplex_scaling = max(simplex_scaling-0.0001f, 1.001f/1.025f );
    }
    if (simplex_num_iters < 5) {
        simplex_scaling = min(simplex_scaling+0.0001f, 1.10f);
    }

    simplex_num_iters = 0;
    // set all negative coordinates to zero.
    to_nonnegative_orthant(x);
    
    // check if the vector is magically already on the simplex
    float s = eval_simplex_disagreement(x, 0, 0);
    //float s = accumulate( x.begin(), x.end(), 0.0f );

    if (abs(s-eta)<tol) { return; }
    // lower and upper bounds are initialized sans-checks.
    // LB=0 is clearly sensible. UB=sum(x) is also, but can be too conservative.
    // ideally one should memorize the solution from the previous iteration,
    // and try to estimate lower/upper bounds based on that; that may save a
    // few iterations.
    
    float scale = base_scale*simplex_scaling;
    float LB = simplex_proj_sol/scale;
    float UB = simplex_proj_sol*scale;
    float val = 0;
    int itr = 0;
    float fval_x, fval_ub, fval_lb;
    float LBog = LB;
    float UBog = UB; 
    int maxiter = 1000;
    fval_ub = eval_simplex_disagreement(x, UB, eta);
    fval_lb = eval_simplex_disagreement(x, LB, eta);
    while (fval_lb < 0) {
        simplex_scaling += 0.0001f;
        LB -= s/500;
        fval_lb = eval_simplex_disagreement(x, LB, eta);
    }
    while (fval_ub > 0) { 
        UB += s/500;
        simplex_scaling += 0.0001f;
        fval_ub = eval_simplex_disagreement(x, UB, eta);
    }
    
    while ( (itr < maxiter) && ((UB-LB)>tol) ) {
        val = (UB+LB)/2.0f;
        fval_x = eval_simplex_disagreement(x, val, eta);
        if (fval_x*fval_ub > 0) {
            UB = val;
            fval_ub = fval_x;
            continue;
        }
        if (fval_x*fval_lb > 0) {
            LB = val;
            fval_lb = fval_x;
            continue;
        }
        if (abs(fval_x) < tol) { break; }
        itr++;
    }
    if (itr >= maxiter) {
        printf("ERROR.\n");
        printf("project_onto_prob_simplex could not find a solution!\n");
        printf("started with lb,ub = [%f,%f], vals = [%f,%f]", LBog, UBog,
                eval_simplex_disagreement(x, LBog, eta),
                eval_simplex_disagreement(x, UBog, eta) );
        printf("Current estimated simplex_disagreement: %f\n", fval_x);
        printf("Will apply this as the solution.\n");
    }
    
    #ifdef HAVE_SSE
    float* x_ = (float*)&x[0];
    int nloops = x.size()/4;
    __m128 m_x;
    __m128 m_zero = _mm_set1_ps(0.0f);
    __m128 m_alpha = _mm_set1_ps(val);
    for (int i = 0; i < nloops; ++i) {
        m_x = _mm_max_ps( _mm_sub_ps( _mm_loadu_ps(x_+4*i), m_alpha ), m_zero );
        _mm_storeu_ps(x_+4*i, m_x);
    }
    for (size_t i = 4*nloops; i < x.size(); ++i) { 
        x[i] = max(x[i]-val, 0.0f);
    }
    #else
    for (size_t i = 0; i < x.size(); ++i) { x[i] = max(0.0f, x[i] - val); }
    #endif
    simplex_proj_sol = val;
    
}

float solver_primaldual::eval_simplex_disagreement(
        const vector<float>& x, const vector<float>& alpha, float eta) {
    simplex_num_iters++;
    #ifdef HAVE_SSE
    return eval_simplex_disagreement_sse(x, alpha, eta);
    #else
    float r = -eta;
    for (size_t i = 0; i < x.size(); ++i) { r += max(x[i]-alpha[i], 0.0f); }
    return r;
    #endif
}

float solver_primaldual::eval_simplex_disagreement(
        const vector<float>& x, float alpha, float eta) {
    simplex_num_iters++;
    #ifdef HAVE_SSE
    return eval_simplex_disagreement_sse(x, alpha, eta);
    #else
    float r = -eta;
    for (size_t i = 0; i < x.size(); ++i) { r += max(x[i]-alpha, 0.0f); }
    return r;
    #endif
}

#ifdef HAVE_SSE
float solver_primaldual::eval_simplex_disagreement_sse(
        const vector<float>& x, float alpha, float eta) { 
    float r = -eta;
    //
    float* x_ = (float*)&x[0];
    int nloops = x.size()/4;
    __m128 mx;
    __m128 m_zero = _mm_set1_ps(0.0f);
    __m128 m_alpha = _mm_set1_ps(alpha);
    __m128 m_sum = _mm_set1_ps(0.0f);
    
    for (int i = 0; i < nloops; ++i) {
        mx = _mm_loadu_ps(x_ + 4*i);
        mx = _mm_sub_ps(mx, m_alpha);
        mx = _mm_max_ps(mx, m_zero);
        m_sum = _mm_add_ps(m_sum, mx);
    }
    float sum[4];
    _mm_storeu_ps(sum, m_sum);
    r += (sum[0] + sum[1] + sum[2] + sum[3]);
    for (size_t i = 4*nloops; i < x.size(); ++i) { 
        r += max(x[i]-alpha, 0.0f);
    }
    return r;
}
#endif


#ifdef HAVE_SSE
float solver_primaldual::eval_simplex_disagreement_sse(
        const vector<float>& x, const vector<float>& alpha, float eta) { 
    float r = -eta;
    //
    float* x_ = (float*)&x[0];
    float* alpha_ = (float*)&alpha[0];
    int nloops = x.size()/4;
    __m128 mx, m_alpha;
    __m128 m_zero = _mm_set1_ps(0.0f);
    __m128 m_sum = _mm_set1_ps(0.0f);
    
    for (int i = 0; i < nloops; ++i) {
        m_alpha = _mm_loadu_ps(alpha_ + 4*i);
        mx = _mm_loadu_ps(x_ + 4*i);
        mx = _mm_sub_ps(mx, m_alpha);
        mx = _mm_max_ps(mx, m_zero);
        m_sum = _mm_add_ps(m_sum, mx);
    }
    float sum[4];
    _mm_storeu_ps(sum, m_sum);
    r += (sum[0] + sum[1] + sum[2] + sum[3]);
    for (size_t i = 4*nloops; i < x.size(); ++i) { 
        r += max(x[i]-alpha[i], 0.0f);
    }
    return r;
}
#endif


void solver_primaldual::solve_y1(vector<float>& DX, float sigma_y) { 
#ifdef HAVE_SSE
    solve_y1_sse(DX, sigma_y);
#else
    for (size_t i = 0; i < sol_y1.size(); ++i) { 
        float y1_hat = sol_y1[i] + sigma_y*DX[i];
        int s = y1_hat>0?-1:+1;
        sol_y1[i] = s*max( -params->prob.weights[i], s*y1_hat );
    }
#endif
}

#ifdef HAVE_SSE
void solver_primaldual::solve_y1_sse(vector<float>& DX, float sigma_y) {
    size_t n = sol_y1.size();
    size_t nloops = n/4;
        
    float* sol_y_ = (float*)&sol_y1[0];
    float* DX_ = (float*)&DX[0];
    float* weights = (float*)&params->prob.weights[0];

    __m128 dx, y, yhat, w, w_neg, temp, t1, t2;

    __m128 m_sigma_y = _mm_set1_ps(sigma_y);
    __m128 m_zero = _mm_set1_ps(0.0f);

    for (size_t i = 0; i < nloops; ++i) {
        dx = _mm_loadu_ps(DX_ + 4*i);
        y = _mm_loadu_ps(sol_y_ + 4*i);
        w = _mm_loadu_ps(weights + 4*i);
        w_neg = _mm_sub_ps(m_zero, w);

        temp = _mm_mul_ps(m_sigma_y, dx); // sigma*DX[i]
        yhat = _mm_add_ps(y, temp); // yhat = y[i] + sigma*DX[i]
        
        // get absolute value in a very convoluted way:
        temp = _mm_min_ps(m_zero, yhat);
        t1 = _mm_max_ps(w_neg, temp);
        temp = _mm_max_ps(m_zero, yhat);
        t2 = _mm_min_ps(w, temp);
        yhat = _mm_add_ps(t1, t2);

        _mm_storeu_ps(sol_y_ + 4*i, yhat);
    }
    for (size_t i = 4*nloops; i < n; ++i) {
        float y = sol_y1[i] + sigma_y*DX[i];
        float w = params->prob.weights[i];
        sol_y1[i] = max(-w, min(0.0f, y)) + min(+w, max(0.0f, y));
    }
}
#endif

void solver_primaldual::solve_y3(float sigma) {
    #ifdef HAVE_SSE
    solve_y3_sse(sigma);
    #else
    for (size_t i = 0; i < sol_y3.size(); ++i) { 
        sol_y3[i] = sol_y3[i] + sigma * sol_x[i];
        
        sol_y3[i] = min( max( min(0.0f, sol_y3[i]-sigma), 
                    -params->prob.kappa[i]), sol_y3[i] );
    }
    #endif
}

#ifdef HAVE_SSE
void solver_primaldual::solve_y3_sse(float sigma) {
    int n = sol_y3.size();
    int nloops = n / 4;
   
    float* sol_y_ = (float*)&sol_y3[0];
    float* sol_x_ = (float*)&sol_x[0];
    float* kappa_ = (float*)&params->prob.kappa[0];

    __m128 m_temp, m_y, m_x, m_kappa, m_negkappa;
    __m128 m_sigma = _mm_set1_ps(sigma);
    __m128 m_zero = _mm_set1_ps(0.0f);

    for (int i = 0; i < nloops; ++i) {
        m_y = _mm_loadu_ps(sol_y_ + 4*i);
        m_x = _mm_loadu_ps(sol_x_ + 4*i);
        m_kappa = _mm_loadu_ps(kappa_ + 4*i);
        m_negkappa = _mm_sub_ps(m_zero, m_kappa);

        m_y = _mm_add_ps( m_y, _mm_mul_ps(m_x, m_sigma));
        m_temp = _mm_sub_ps(m_y, m_sigma);
        m_temp = _mm_min_ps( m_y, 
                _mm_max_ps( _mm_min_ps(m_zero, m_temp), m_negkappa) );
        _mm_storeu_ps( sol_y_ + 4*i, m_temp );
    }
    
    for (size_t i = 4*nloops; i < sol_y3.size(); ++i) { 
        sol_y3[i] = sol_y3[i] + sigma * sol_x[i];
        sol_y3[i] = min( sol_y3[i], 
                max( sol_y3[i] - sigma, -params->prob.kappa[i] ));
    }
}
#endif

float solver_primaldual::solve_c_ell_infinity(
        vector<float>& DTY, vector<float>& DOCCTY, float sigma_c) {
    vector<float> c_temp = vector<float>( sol_c.size() );
    #ifdef HAVE_SSE
        solve_c_ell_infinity_prepare_(c_temp, DTY, DOCCTY, sigma_c, params->opt.use_temporal_penalty );
    #else
        float ub = 0;
        if (params->opt.layer_upper_bound > 0) { ub = params->opt.layer_upper_bound; }
        if (params->opt.layer_upper_bound <= 0) { ub = 9999; }

        if (params->opt.use_temporal_penalty) { 
            // do something
            for (size_t i = 0; i < sol_c.size(); ++i) { 
                c_temp[i] = max(0.0f, sol_c[i] - sigma_c*(DTY[i] + DOCCTY[i] + sol_y3[i]) );
            }
        } else {
            // this can easily be speeded up significantly.....
            for (size_t i = 0; i < sol_c.size(); ++i) {
                c_temp[i] = max(0.0f, sol_c[i] - sigma_c*(DTY[i] + DOCCTY[i]));
            }
        }
    #endif

    // if sum(max(0,x))>1, then need to project onto a simplex.
    float s = eval_simplex_disagreement(c_temp, 0, 0);
    vector<float> c_simplex = c_temp;
    if (s > sigma_c*params->prob.tau2) {
        project_onto_prob_simplex(c_simplex, sigma_c*params->prob.tau2);
    }
    // update the solution:
    float dx_ = 0.0f;
    #ifdef HAVE_SSE
    dx_ = solve_c_ell_infinity_update_(c_temp, c_simplex);
    #else
    float c_prev, dc;
    for (size_t i = 0; i < sol_c.size(); ++i) { 
        c_prev = sol_c[i];
        sol_c[i] = max(0.0f, c_temp[i] - c_simplex[i]);
        dc = sol_c[i]-c_prev;
        sol_x[i] = sol_c[i] + dc;
        dx_ += dc*dc;
    }
    #endif
    return dx_;
}

float solver_primaldual::solve_c_ell1_ell_infinity(
        vector<float>& DTY, vector<float>& DOCCTY, float sigma_c) { 
    float sigma_inv = 1.0f/sigma_c;
    size_t nloops = sol_c.size()/4;

    // FIRST PREPARE THE VARIABLE (NOTE HERE IS NO PROJECTION INTO FEASIBLE SET) 
    vector<float> c_hat = vector<float>( sol_c.size() );
     
    #ifdef HAVE_SSE
    int project_onto_nonnegative_orthant = 0;
    solve_c_ell_infinity_prepare_(c_hat, DTY, DOCCTY, sigma_c, 
            params->opt.use_temporal_penalty, project_onto_nonnegative_orthant);
    #else
    if (params->opt.use_temporal_penalty) { 
        for (size_t i = 0; i < sol_c.size(); ++i) { 
            c_hat[i] = sol_c[i] - sigma_c*(DTY[i] + DOCCTY[i] + sol_y3[i]);
        }
    } else {
        for (size_t i = 0; i < sol_c.size(); ++i) { 
            c_hat[i] = sol_c[i] - sigma_c*(DTY[i] + DOCCTY[i]);
        }
    }
    #endif
    
    vector<float> c_temp = c_hat;
    vector<float> c_proj;
    // PREMULTIPLY BY INVERSE SIGMA
    #ifdef HAVE_SSE
    __m128 m_c;
    __m128 m_sigmainv = _mm_set1_ps(sigma_inv);

    float* c_temp_ = &c_temp[0];
    for (size_t i = 0; i < nloops; ++i) { 
        m_c = _mm_loadu_ps(c_temp_+4*i);
        _mm_storeu_ps(c_temp_+4*i, _mm_mul_ps(m_c, m_sigmainv) );
    }
    for (size_t i = 4*nloops; i < c_temp.size(); ++i) { c_temp[i]*= sigma_inv; }
    #else
    for (size_t i = 0; i < c_temp.size(); ++i) { c_temp[i]*= sigma_inv; }
    #endif
    
    // NOW PROJECT:
    float s = eval_simplex_disagreement( c_temp, params->prob.tau1, 0.0f );
    bool in_simplex = (s <= params->prob.tau2);
    if (in_simplex) { 
        c_proj = vector<float>( sol_c.size() );
        for (size_t i = 0; i < c_proj.size(); ++i) { 
            c_proj[i] = c_hat[i];            
        }
        for (size_t i = 0; i < c_temp.size(); ++i) { 
            sol_c[i] = 0; 
        }
    } else {
        #ifdef HAVE_SSE
        __m128 m_c,m_tau1;
        float* c_temp_ = &c_temp[0];
        float* c_tau = &params->prob.tau1[0];

        for (size_t i = 0; i < nloops; ++i) { 
            m_tau1 = _mm_loadu_ps(c_tau + 4*i);
            m_c = _mm_loadu_ps(c_temp_+4*i);
            _mm_storeu_ps( c_temp_+4*i, _mm_sub_ps(m_c, m_tau1));
        }
        for (size_t i = 4*nloops; i < c_temp.size(); ++i) {
            c_temp[i] -= params->prob.tau1[i];            
        }
        #else

        for (size_t i = 0; i < c_temp.size(); ++i) { c_temp[i] = c_temp[i] - params->prob.tau1[i]; }
        #endif
        // project this onto the simplex:
        c_proj = c_temp;
        project_onto_prob_simplex(c_proj, params->prob.tau2 );
        
        // here it is a bit difficult to put it into SSE...
        // recover negative entries
        ell1_ell_infinity_recover_neg(c_proj, c_temp, sigma_c, params->prob.tau1);
    }
    #ifdef HAVE_SSE
    return solve_c_ell_infinity_update_(c_hat, c_proj);
    #else
    float ub = 0;
    if (params->opt.layer_upper_bound > 0) { ub = params->opt.layer_upper_bound; }
    if (params->opt.layer_upper_bound <= 0) { ub = 9999; }
    float dx = 0;
    float c_prev, dc;
    for (size_t i = 0; i < sol_c.size(); ++i) { 
        c_prev = sol_c[i];
        sol_c[i] = min(ub, max(0.0f, c_hat[i] - c_proj[i]));
        dc = sol_c[i]-c_prev;
        sol_x[i] = sol_c[i] + dc;
        dx += dc*dc;
    }
    return dx;
    #endif
}

void solver_primaldual::ell1_ell_infinity_recover_neg(vector<float>& c_out, 
                                const vector<float>& c, float sigma, const vector<float>& tau) {
    #ifdef HAVE_SSE
    //----------------------------------------------------
    float* c_out_ = &c_out[0];
    const float* c_ = &c[0];
    const float* tau_ = &tau[0];
    __m128 m_c, m_cout, m_tau, mask, notmask;
    
    __m128 m_one = _mm_set1_ps(1.0f);    
    __m128 m_zero = _mm_set1_ps(0.0f);
    __m128 m_sigma = _mm_set1_ps(sigma);
    __m128 m_ffff = _mm_cmpgt_ps( m_one, m_zero );

    size_t nloops = c.size()/4;

    for (size_t i = 0; i < nloops; ++i) {
        m_c = _mm_loadu_ps( c_ + 4*i);
        m_cout = _mm_loadu_ps( c_out_ + 4*i);
        m_tau = _mm_loadu_ps( tau_ + 4*i);

        mask = _mm_cmplt_ps(m_c, m_zero); // c(i)<0
        notmask = _mm_xor_ps( mask, m_ffff ); // c(i)>= 0
        m_cout = _mm_or_ps( _mm_and_ps( mask, m_c), _mm_and_ps( notmask, m_cout) );

        m_cout = _mm_mul_ps( _mm_add_ps( m_cout, m_tau ), m_sigma );
        _mm_storeu_ps(c_out_ + 4*i, m_cout);
    }
    for (size_t i = nloops*4; i < c.size(); ++i) {
        if (c[i] < 0) { c_out[i] = c[i]; }
        c_out[i] = c_out[i] + tau[i]; 
        c_out[i] *= sigma;
    }
    //----------------------------------------------------
    #else
        for (size_t i = 0; i < c.size(); ++i) {
            if (c[i] < 0) { c_out[i] = c[i]; }
            c_out[i] = c_out[i] + tau[i]; 
            c_out[i] *= sigma;
        }
    #endif    
}



#ifdef HAVE_SSE
void solver_primaldual::solve_c_ell_infinity_prepare_(vector<float>& out, 
        vector<float>& DTY, vector<float>& DOCCTY, float sigma, bool temporal,
        int project_onto_feasible_set) { 
    float* out_ = (float*)&out[0];        
    float* DTY_ = (float*)&DTY[0];        
    float* DOCCTY_ = (float*)&DOCCTY[0];        
    float* sol_y3_ = (float*)&sol_y3[0];
    float* sol_c_ = (float*)&sol_c[0];

    int n = sol_c.size();
    int nloops = n / 4;
    
    __m128 m_sigma = _mm_set1_ps(sigma);
    __m128 m_temp, m_c, m_dty, m_doccty, m_y3;
    if (temporal == 0) { // without temporal penalty (i.e. without y3)
        for (int i = 0; i < nloops; ++i) { 
            m_c = _mm_loadu_ps( sol_c_ + 4*i);
            m_dty = _mm_loadu_ps( DTY_ + 4*i);
            m_doccty = _mm_loadu_ps( DOCCTY_ + 4*i);

            m_temp = _mm_mul_ps( m_sigma, _mm_add_ps(m_dty, m_doccty) );          
            m_temp = _mm_sub_ps( m_c, m_temp );

            _mm_storeu_ps(out_ + 4*i, m_temp);
        }
        for (int i = nloops*4; i < n; ++i) { 
            out[i] = sol_c[i] - sigma*(DTY[i]+DOCCTY[i]);
        }
    } else { // with temporal penalty (i.e. with y3)
        for (int i = 0; i < nloops; ++i) { 
            m_c = _mm_loadu_ps( sol_c_ + 4*i);
            m_dty = _mm_loadu_ps( DTY_ + 4*i);
            m_doccty = _mm_loadu_ps( DOCCTY_ + 4*i);
            m_y3 = _mm_loadu_ps( sol_y3_ + 4*i);

            m_temp = _mm_mul_ps( m_sigma, _mm_add_ps(m_y3, _mm_add_ps(m_dty, m_doccty)) );          
            m_temp = _mm_sub_ps( m_c, m_temp );

            _mm_storeu_ps(out_ + 4*i, m_temp);
        }
        for (int i = nloops*4; i < n; ++i) { 
            out[i] = sol_c[i] - sigma*(DTY[i]+DOCCTY[i]+sol_y3[i]);
        }
    }
    if (project_onto_feasible_set) { 
       to_nonnegative_orthant(out); 
    }
}

float solver_primaldual::solve_c_ell_infinity_update_(
        vector<float>& c_temp, vector<float>& c_simplex) { 
    float* c_simplex_ = (float*)&c_simplex[0];        
    float* c_temp_ = (float*)&c_temp[0];
    float* sol_c_ = (float*)&sol_c[0];
    float* sol_x_ = (float*)&sol_x[0];

    float ub = 0;
    if (params->opt.layer_upper_bound > 0) { ub = params->opt.layer_upper_bound; }
    if (params->opt.layer_upper_bound <= 0) { ub = 9999; }

    int n = sol_c.size();
    int nloops = n / 4;
    
    __m128 m_zero = _mm_set1_ps(0.0f);
    __m128 m_ub = _mm_set1_ps(ub);

    __m128 m_temp, m_c, m_c_temp, m_cprev, m_dc, m_x, m_c_simplex;
    __m128 m_dx = _mm_set1_ps(0.0f);

    for (int i = 0; i < nloops; ++i) { 
        m_c = _mm_loadu_ps( sol_c_ + 4*i);
        m_c_simplex = _mm_loadu_ps( c_simplex_ + 4*i);
        m_c_temp = _mm_loadu_ps( c_temp_ + 4*i);
        //
        m_cprev = m_c;
        m_temp = _mm_sub_ps( m_c_temp , m_c_simplex );
        m_c = _mm_min_ps( _mm_max_ps(m_zero, m_temp), m_ub );
        m_dc = _mm_sub_ps( m_c, m_cprev );
        m_x = _mm_add_ps( m_c, m_dc );
        m_dx = _mm_add_ps( m_dx, _mm_mul_ps(m_dc, m_dc ) );    

        _mm_storeu_ps(sol_x_ + 4*i, m_x);
        _mm_storeu_ps(sol_c_ + 4*i, m_c);
    }
    // fix dx:  
    float dxsum_[4];
    _mm_storeu_ps(dxsum_, m_dx);
    float dx_ = dxsum_[0] + dxsum_[1] + dxsum_[2] + dxsum_[3];
     
    for (int i = nloops*4; i < n; ++i) { 
        for (size_t i = 0; i < sol_c.size(); ++i) { 
            float c_prev = sol_c[i];
            sol_c[i] = min(ub, max(0.0f, c_temp[i] - c_simplex[i]));
            float dc = sol_c[i]-c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
        }
    }
    return dx_;
}
#endif


float solver_primaldual::solve_c(vector<float>& DTY, 
        vector<float>& DOCCTY, float sigma, vector<float>& tau) { 
    
    #ifdef HAVE_SSE
    return solve_c_sse(DTY, DOCCTY, sigma, tau);
    #else
    // upper bound. this is a little bit faster.
    float c_prev;
    float dx_ = 0.0f;
    float dc;
    if (params->opt.layer_upper_bound > 0) { 
        for (size_t i = 0; i < sol_c.size(); ++i) {
            c_prev = sol_c[i];
            sol_c[i] = min(params->opt.layer_upper_bound,
                    max(0.0f, sol_c[i] - sigma*(DTY[i] + DOCCTY[i] + tau[i]) ) );
            #ifdef NO_DX_CHECK
            sol_x[i] = 2*sol_c[i] - c_prev;
            #else
            dc = sol_c[i] - c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
            #endif
        }
    } else {
        for (size_t i = 0; i < sol_c.size(); ++i) {
            c_prev = sol_c[i];
            sol_c[i] = max(0.0f, sol_c[i] - sigma*(DTY[i] + DOCCTY[i] + tau[i]) );
            #ifdef NO_DX_CHECK
            sol_x[i] = 2*sol_c[i] - c_prev;
            #else
            dc = sol_c[i] - c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
            #endif
        }
    }
    #ifdef NO_DX_CHECK
    return 1.0f;
    #else
    return dx_;
    #endif
    #endif
}

float solver_primaldual::solve_c_temporal(vector<float>& DTY, 
        vector<float>& DOCCTY, float sigma, vector<float>& tau) { 
    
    #ifdef HAVE_SSE
    return solve_c_temporal_sse(DTY, DOCCTY, sigma, tau);
    #else
    // use a separate block for the case that hte user specified an
    // upper bound. this is a little bit faster.
    float c_prev;
    float dx_ = 0.0f;
    float dc;

    if (params->opt.layer_upper_bound > 0) {
        for (size_t i = 0; i < sol_c.size(); ++i) {
            float sigmatau = sigma*tau[i];
            float oneplussigmatau = 1+sigma*tau[i];
            c_prev = sol_c[i];
            // ----------------------------------------------------------------
            float z = sol_c[i] - sigma*(DTY[i]+DOCCTY[i]);
            if (z > oneplussigmatau) {
                sol_c[i] = min(z - sigmatau, params->opt.layer_upper_bound);
            } else if (z < oneplussigmatau - sigma*params->prob.kappa[i]) { 
                sol_c[i] = max(0.0f, z - sigmatau + sigma*params->prob.kappa[i]);
            } else {
                sol_c[i] = 1;
            }
            // ----------------------------------------------------------------
            #ifdef NO_DX_CHECK
            sol_x[i] = 2*sol_c[i] - c_prev;
            #else
            dc = sol_c[i] - c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
            #endif
        }
    } else {
        for (size_t i = 0; i < sol_c.size(); ++i) {
            float sigmatau = sigma*tau[i];
            float oneplussigmatau = 1+sigma*tau[i];

            c_prev = sol_c[i];
            // ----------------------------------------------------------------
            float z = sol_c[i] - sigma*(DTY[i]+DOCCTY[i]);
            if (z > oneplussigmatau) {
                sol_c[i] = z - sigmatau;
            } else if (z < oneplussigmatau - sigma*params->prob.kappa[i]) { 
                sol_c[i] = max(0.0f, z - sigmatau + sigma*params->prob.kappa[i]);
            } else {
                sol_c[i] = 1;
            }
            // ----------------------------------------------------------------
            #ifdef NO_DX_CHECK
            sol_x[i] = 2*sol_c[i] - c_prev;
            #else
            dc = sol_c[i] - c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
            #endif
        }
    }
    #ifdef NO_DX_CHECK
    dx_ = 1.0f;
    #endif
    return dx_;
    #endif
}


#ifdef HAVE_SSE
float solver_primaldual::solve_c_sse(vector<float>& DTY,
        vector<float>& DOCCTY, float sigma, vector<float>& tau_) {
    size_t n = sol_c.size();
    size_t nloops = n/4;
    __m128 c, c_prev, x, dty, doccty, dc, dc2, temp, m_tau;
    __m128 dxsum = _mm_set1_ps(0.0f);

    __m128 m_sigma = _mm_set1_ps(sigma);
    __m128 m_zero = _mm_set1_ps(0.0f);

    float* DTY_ = &DTY[0];
    float* DOCCTY_ = &DOCCTY[0];
    float* sol_c_ = &sol_c[0];
    float* sol_x_ = &sol_x[0];
    float* tau = &tau_[0];

    if (params->opt.layer_upper_bound > 0) {
        __m128 m_ub = _mm_set1_ps(params->opt.layer_upper_bound);
        for (size_t i = 0; i < nloops; i++) { 
            m_tau = _mm_loadu_ps(tau + 4*i);
            c = _mm_loadu_ps(sol_c_ + 4*i);
            dty = _mm_loadu_ps(DTY_ + 4*i);
            doccty = _mm_loadu_ps(DOCCTY_ + 4*i);
            c_prev = c;

            temp = _mm_add_ps(dty, doccty);
            temp = _mm_add_ps(temp, m_tau);
            temp = _mm_mul_ps(temp, m_sigma); // sigma_c*(DTY+DOCCTY+tau)
            c = _mm_sub_ps(c, temp);
            c  = _mm_max_ps(c, m_zero); // max(0, c)
            c  = _mm_min_ps(c, m_ub);
            dc = _mm_sub_ps(c, c_prev);

            x = _mm_add_ps(c, dc); // x = sol_c + dc
            _mm_storeu_ps(sol_c_ + 4*i, c);
            _mm_storeu_ps(sol_x_ + 4*i, x);
           
            #ifndef NO_DX_CHECK 
            dc2 = _mm_mul_ps(dc, dc);
            dxsum = _mm_add_ps(dxsum, dc2);
            #endif
        }
    } else {
        for (size_t i = 0; i < nloops; i++) { 
            m_tau = _mm_loadu_ps(tau + 4*i);            
            c = _mm_loadu_ps(sol_c_ + 4*i);
            dty = _mm_loadu_ps(DTY_ + 4*i);
            doccty = _mm_loadu_ps(DOCCTY_ + 4*i);
            c_prev = c;

            temp = _mm_add_ps(dty, doccty);
            temp = _mm_add_ps(temp, m_tau);
            temp = _mm_mul_ps(temp, m_sigma); // sigma_c*(DTY+DOCCTY+tau)
            c = _mm_sub_ps(c, temp);
            c  = _mm_max_ps(c, m_zero); // max(0, c)

            dc = _mm_sub_ps(c, c_prev);
            x = _mm_add_ps(c, dc); // x = sol_c + dc

            _mm_storeu_ps(sol_c_ + 4*i, c);
            _mm_storeu_ps(sol_x_ + 4*i, x);
           
            #ifndef NO_DX_CHECK 
            dc2 = _mm_mul_ps(dc, dc);
            dxsum = _mm_add_ps(dxsum, dc2);
            #endif
        }
    }
    // need to verify:
    float dxsum_[4];
    _mm_storeu_ps(dxsum_, dxsum);
    float dx = dxsum_[0] + dxsum_[1] + dxsum_[2] + dxsum_[3];
    for (size_t i = 4*nloops; i < sol_c.size(); ++i) {
        float c_prev = sol_c[i];
        sol_c[i] = max(0.0f, sol_c[i] - sigma*(DTY[i] + DOCCTY[i] + tau[i]) );
        if (params->opt.layer_upper_bound > 0) {
            sol_c[i] = min(sol_c[i], params->opt.layer_upper_bound);
        }
        float dc = sol_c[i] - c_prev;
        sol_x[i] = sol_c[i] + dc;
        #ifndef NO_DX_CHECK
        dx += dc*dc;
        #endif
    }
    #ifdef NO_DX_CHECK
    dx = 1.0f;
    #endif
    return dx;
}


float solver_primaldual::solve_c_temporal_sse(vector<float>& DTY,
        vector<float>& DOCCTY, float sigma, vector<float>& tau_) {
    size_t n = sol_c.size();
    size_t nloops = n/4;
    
    float dc;

    __m128 c, m_cprev, dty, doccty, m_dc, m_dc2, temp;
    __m128 m_tau, m_sigmatau, m_oneplussigmatau;
    __m128 dxsum = _mm_set1_ps(0.0f);
    __m128 m_sigma = _mm_set1_ps(sigma);
    __m128 m_zero = _mm_set1_ps(0.0f);
    __m128 m_one = _mm_set1_ps(1.0f);
    __m128 m_ffff = _mm_cmpgt_ps( m_one, m_zero);
    float* DTY_ = &DTY[0];
    float* DOCCTY_ = &DOCCTY[0];
    float* sol_c_ = &sol_c[0];
    float* sol_x_ = &sol_x[0];
    float* kappa_ = &params->prob.kappa[0];
    float* tau = &tau_[0];

    __m128 m_ub = _mm_set1_ps( params->opt.layer_upper_bound );
    if (params->opt.layer_upper_bound <= 0) {
        m_ub = _mm_set1_ps( 99999999.0f ); // idgaf, less code.
    }
    
    __m128 m_kappa, m1, m2, m3, out1, out2, out3;

    for (size_t i = 0; i < nloops; ++i) { 
        // ----------------------------------------------------------------
        // first just compute sol_c <- sol_c - sigma*(DTY[i]+DOCCTY[i])
        m_tau = _mm_loadu_ps(tau + 4*i);
        m_sigmatau = _mm_mul_ps(m_sigma, m_tau);
        m_oneplussigmatau = _mm_add_ps(m_one, m_sigmatau);

        c = _mm_loadu_ps(sol_c_ + 4*i);
        m_cprev = c;
        dty = _mm_loadu_ps(DTY_ + 4*i);
        doccty = _mm_loadu_ps(DOCCTY_ + 4*i);
        temp = _mm_add_ps(dty, doccty);
        temp = _mm_mul_ps(temp, m_sigma); // sigma_c*(DTY+DOCCTY)
        c = _mm_sub_ps(c, temp);
        // ----------------------------------------------------------------

        m_kappa = _mm_loadu_ps(kappa_ + 4*i);
        __m128 m_sigmakappa = _mm_mul_ps(m_sigma, m_kappa);
        __m128 m_lobound = _mm_sub_ps(m_oneplussigmatau, m_sigmakappa);

        // m1 = ffff , when c >= 1+sigma*tau
        m1 = _mm_cmpgt_ps(c, m_oneplussigmatau);
        // m2 = ffff , when c <= 1+sigma*tau - sigma*kappa
        m2 = _mm_cmplt_ps(c, m_lobound);
        // m3 = ffff , when m1=0 & m2 = 0
        m3 = _mm_xor_ps( _mm_or_ps(m1, m2), m_ffff );

        // now compute all outputs:
        out1 = _mm_sub_ps(c, m_sigmatau);
        out2 = _mm_add_ps(out1, m_sigmakappa);
        out1 = _mm_min_ps( out1, m_ub );
        out2 = _mm_max_ps( out2, m_zero );
        out3 = m_one;
        // and all three, then or them:
        out1 = _mm_and_ps(out1, m1);
        out2 = _mm_and_ps(out2, m2);
        out3 = _mm_and_ps(out3, m3);
        c = _mm_or_ps( _mm_or_ps(out1, out2), out3);

        _mm_storeu_ps( sol_c_ + 4*i, c);
        //---------------------------------------------------------------------
        m_dc = _mm_sub_ps(c, m_cprev);
        temp = _mm_add_ps(c, m_dc );
        #ifndef NO_DX_CHECK
        m_dc2 = _mm_mul_ps(m_dc, m_dc);
        dxsum = _mm_add_ps(m_dc2, dxsum);
        #endif 
        _mm_storeu_ps( sol_x_ + 4*i, temp);
    }
    float dxsum_[4];
    _mm_storeu_ps(dxsum_, dxsum);
    float dx = dxsum_[0] + dxsum_[1] + dxsum_[2] + dxsum_[3];
     
    for (size_t i = 4*nloops; i < sol_c.size(); ++i) { 
        float sol_cprev = sol_c[i];
        float sigmatau = sigma*tau[i];
        float oneplussigmatau = 1 + sigmatau;
        sol_c[i] = sol_c[i] - sigma*(DTY[i] + DOCCTY[i] );

        if (sol_c[i] > oneplussigmatau) {
            sol_c[i] = min(sol_c[i] - sigmatau, params->opt.layer_upper_bound);
        } else if (sol_c[i] < oneplussigmatau - sigma*params->prob.kappa[i]) { 
            sol_c[i] = max(0.0f, sol_c[i] - sigmatau + sigma*params->prob.kappa[i]);
        } else {
            sol_c[i] = 1;
        }
        #ifdef NO_DX_CHECK
        sol_x[i] = 2*sol_c[i]-sol_cprev;
        #else
        dc = sol_c[i] - sol_cprev;
        sol_x[i] = sol_c[i] + dc;
        dx += dc*dc;
        #endif
    }
    #ifdef NO_DX_CHECK
    dx = 1.0f;
    #endif
    return dx;
}
#endif


void solver_primaldual::to_nonnegative_orthant(vector<float>& x) { 
    size_t n = x.size();
#ifdef HAVE_SSE
    size_t nloops = n / 4;
    __m128 m_x;
    __m128 m_zero = _mm_set1_ps(0.0f);
    float* x_ = &x[0];
    for (size_t i = 0; i < nloops; ++i) { 
        // ----------------------------------------------------------------
        m_x = _mm_max_ps( _mm_loadu_ps(x_ + 4*i), m_zero);
        _mm_storeu_ps(x_ + 4*i, m_x);
    } 
    for (size_t i = 4*nloops; i < x.size(); ++i) { x[i] = max( x[i], 0.0f ); }
#else
    for (size_t i = 0; i < x.size(); ++i) { x[i] = max( x[i], 0.0f ); }
#endif
}

//! approximates \|D^T D\|_2 using the power-method
float solver_primaldual::estimate_divgrad_norm() {
    if (params->prob.solve_pixelwise) { 
        printf("for a pixel lattice, |D^T D|_2 = 8 - no need to estimate that\n");
    }
    int niters = 20;    
    vector<float> x = vector<float>( params->prob.nnodes );
    vector<float> DX = vector<float>( params->prob.nedges );
    for (size_t i = 0; i < x.size(); ++i) { 
        x[i] = (float)(double( rand() ) / double(RAND_MAX)) - 0.5f; // in [-.5, +.5]
    }
    for (int i = 0; i < niters; ++i) { 
        // normalize the vector to have ||x||_2 = 1
        float norm = 0.0;
        for (size_t u = 0; u < x.size(); ++u) { norm += x[u]*x[u]; }
        norm = sqrt(norm);
        transform(x.begin(), x.end(), x.begin(), bind1st(multiplies<float>(), 1.0f/norm ));
        // apply D^T D
        if (params->prob.solve_pixelwise) { 
            apply_pixelwise_gradient_op(DX, x, params->prob.rows, params->prob.cols);
            apply_pixelwise_gradient_op_transpose(x, DX, params->prob.rows, params->prob.cols);
        } else {
            apply_sparse_op<int>(DX, params->prob.D, x );
            apply_sparse_op<int>(x, params->prob.D, DX, 'T');
        }
    }
    float norm = 0.0;
    for (size_t u = 0; u < x.size(); ++u) { norm += x[u]*x[u]; } norm = sqrt(norm);
    return norm;
}

void solver_primaldual::initialize(vector<float>& layers) {
    if (params == NULL) { 
        // frankly at this stage, should not initialize.
        printf("error; parameters not initialized, but initial point is given\n");
        return; 
    }
    if ((size_t)params->prob.nnodes != layers.size() ) {
        printf("number of nodes does not match initial point dimension.\n");
        printf("this is an error.\n");
        return;
    }
    sol_x = layers;
    sol_c = layers;
}

int solver_primaldual::check_params() {
    if (params == NULL) {
        printf("solver_primaldual: parameters were not set.\n");
        return 0;
    }
    int err = 0;
    err = err | (params->opt.sigma_y < 0);
    err = err | (params->opt.sigma_c < 0);
    if (err) {
        printf("solver_primaldual: optimization step-sizes invalid\n");
        return 0;
    }
    err = err | (params->opt.max_iterations < 1);
    err = err | (params->opt.fx_tolerance < 0); 
    err = err | (params->opt.dx_tolerance < 0); 
    if (err) {
        printf("solver_primaldual: termination conditions invalid.\n");
        return 0;
    }
    err = (params->prob.tau1.size()==0) || (params->prob.tau2 < 0);
    if (err) { 
        printf("solver_primaldual: problem weighing wrong (tau)\n");
        return 0;
    }
    err = err | (params->prob.nnodes <= 0);
    err = err | (params->prob.nedges <= 0);
    err = err | (params->prob.nocc_constraints < 0);
    if (err) { 
        printf("solver_primaldual: problem size set incorrectly\n");
        printf("number of nodes: %d\n", params->prob.nnodes);
        printf("number of edges: %d\n", params->prob.nedges);
        printf("number of constraints: %d\n", params->prob.nocc_constraints);
        return 0;
    }
    if (params->prob.solve_pixelwise == 0) {
        int m1 = params->prob.D.nrows;
        int m2 = params->prob.D.ncols;
        if ((params->prob.nedges != m1) || (params->prob.nnodes != m2)) {
            printf("chosen to use sparse.ops (D) which is not set.\n");
            return 0;
        }
    }
    return 1;
}
