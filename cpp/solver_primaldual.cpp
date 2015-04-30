#include "solver_primaldual.h"

solver_primaldual::solver_primaldual() { 
    params = NULL;  
}

solver_primaldual::~solver_primaldual() {
    if (sol_x)  { free(sol_x); }
    if (sol_c)  { free(sol_c); }
    if (sol_y1) { free(sol_y1); }
    if (sol_y2) { free(sol_y2); }

    if (kappa) { free(kappa); }
    if (tau1) { free(tau1); }
    if (weights) { free(weights); }
    if (occweights) { free(occweights); }
}

solver_primaldual::solver_primaldual(primaldual_params* p) {
    // stub
    params = p;
    sol_x = NULL;
    sol_c = NULL;
    sol_y1 = NULL;
    sol_y2 = NULL;

    kappa = NULL;
    tau1 = NULL;
    weights = NULL;
    occweights = NULL;

    // initialize variables:
    posix_memalign( (void**)&sol_y1, MEM_ALIGN_SIZE, sizeof(float)*params->prob.nedges );
    posix_memalign( (void**)&sol_y2, MEM_ALIGN_SIZE, sizeof(float)*params->prob.nocc_constraints );
    posix_memalign( (void**)&sol_x , MEM_ALIGN_SIZE, sizeof(float)*params->prob.nnodes );
    posix_memalign( (void**)&sol_c , MEM_ALIGN_SIZE, sizeof(float)*params->prob.nnodes );

    if (params->opt.use_temporal_penalty) {
        posix_memalign( (void**)&kappa, MEM_ALIGN_SIZE, sizeof(float)*params->prob.nnodes );
        //kappa = (float*)aligned_alloc( MEM_ALIGN_SIZE, sizeof(float)*params->prob.nnodes );
        memcpy(kappa, &params->prob.kappa[0], sizeof(float)*params->prob.nnodes );   
    }
    posix_memalign( (void**)&tau1,      MEM_ALIGN_SIZE, sizeof(float)*params->prob.nnodes );
    posix_memalign( (void**)&weights,   MEM_ALIGN_SIZE, sizeof(float)*params->prob.nedges );
    posix_memalign( (void**)&occweights,MEM_ALIGN_SIZE, sizeof(float)*params->prob.nocc_constraints );

    memcpy(tau1, &params->prob.tau1[0], sizeof(float)*params->prob.nnodes );
    memcpy(weights, &params->prob.weights[0], sizeof(float)*params->prob.nedges );
    memcpy(occweights, &params->prob.occweights[0], sizeof(float)*params->prob.nocc_constraints );
}

vector<float> solver_primaldual::solve() {
    if (!check_params()) {
        return vector<float>();
    }
    memset(sol_y1,  0, sizeof(float)*params->prob.nedges );
    memset(sol_y2,  0, sizeof(float)*params->prob.nocc_constraints );
    memset(sol_x,   0, sizeof(float)*params->prob.nnodes );
    memset(sol_c,   0, sizeof(float)*params->prob.nnodes );
    if (params->opt.use_initial_layers) { 
        memcpy( sol_c, &params->opt.init_layers[0], sizeof(float)*params->prob.nnodes ); 
        memcpy( sol_x, &params->opt.init_layers[0], sizeof(float)*params->prob.nnodes ); 
    }
       
    // temporary (intermediate) variables
    float *DX, *DOCCX, *DTY, *DOCCTY;
    posix_memalign( (void**)&DX,        MEM_ALIGN_SIZE, sizeof(float)*params->prob.nedges );
    posix_memalign( (void**)&DOCCX,     MEM_ALIGN_SIZE, sizeof(float)*params->prob.nocc_constraints );
    posix_memalign( (void**)&DTY,       MEM_ALIGN_SIZE, sizeof(float)*params->prob.nnodes );
    posix_memalign( (void**)&DOCCTY,    MEM_ALIGN_SIZE, sizeof(float)*params->prob.nnodes );
    
    memset( DX,     0, sizeof(float)*params->prob.nedges );
    memset( DOCCX,  0, sizeof(float)*params->prob.nocc_constraints );
    memset( DTY,    0, sizeof(float)*params->prob.nnodes );
    memset( DOCCTY, 0, sizeof(float)*params->prob.nnodes );
    
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
        for (size_t i = 0; i < params->prob.nocc_constraints; ++i) {
            float y2_hat = sol_y2[i] + sigma_y * DOCCX[i];
            sol_y2[i] = min( max( y2_hat - sigma_y, -params->prob.occweights[i] ), 0.0f);
        }
        if (iter % params->opt.verbosity == 0) { elapsed_y2 = (now_ms()-t2); }
        
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

        // use L1 regularization on the layer values. this is fast.
        if (params->opt.use_temporal_penalty) { 
            dx_ = solve_c_temporal(DTY, DOCCTY, sigma_c, tau1 );
        } else {
            dx_ = solve_c(DTY, DOCCTY, sigma_c, tau1 );
        }

        if (iter % params->opt.verbosity == 0) {
            elapsed_c = (now_ms()-t2);
            elapsed = (now_ms()-t1);
        }
        dx_vec[iter % DFX_SIZE] = sqrt(dx_)/float(params->prob.nnodes);
        dx = accumulate(dx_vec.begin(), dx_vec.end(), 0.0f)/float(DFX_SIZE); 
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
        if ((iter>0) && (params->opt.verbosity>0) && (iter % params->opt.verbosity == 0)) { 
            printf("%4d: fx=%.4e | dfx=%.4e | dx=%.4e | %.3f ms | %.3f ms ||", 
                    iter, fx, dfx, dx, elapsed, fcalc_elapsed);
            printf(" %.3f | %.3f | %.3f |", elapsed_y1, elapsed_y2, elapsed_c);
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

    free(DX);
    free(DOCCX);
    free(DTY);
    free(DOCCTY);
    // ------------------------------------------------------------------------
    vector<float> sol_c_vec = vector<float>(params->prob.nnodes );
    memcpy( &sol_c_vec[0], sol_c, sizeof(float)*params->prob.nnodes );
    // ------------------------------------------------------------------------
    return sol_c_vec;
}

float solver_primaldual::eval_opt_func_value() {
    float *Doccx, *Dx;
    posix_memalign((void**)&Doccx, MEM_ALIGN_SIZE, sizeof(float)*params->prob.nocc_constraints );
    posix_memalign((void**)&Dx,    MEM_ALIGN_SIZE, sizeof(float)*params->prob.nedges );

    if (params->prob.solve_pixelwise) {
        apply_pixelwise_gradient_op(Dx, sol_c, params->prob.rows, params->prob.cols);
    } else {
        apply_sparse_op<int>(Dx, params->prob.D, sol_c);
    }
    apply_sparse_op<int>(Doccx, params->prob.Docc, sol_c);
    //
    float fx = 0.0f;
    // \tau 1^T c
    for (size_t i = 0; i < params->prob.nnodes; ++i) { fx += tau1[i] * sol_c[i]; }
    // \lambda^T max(0, 1 - Docc x)
    for (size_t i = 0; i < params->prob.nocc_constraints; ++i) { 
        fx += occweights[i] * max( 0.0f, 1.0f - Doccx[i] );
    }
    // |WDx|_1
    for (size_t i = 0; i < params->prob.nedges; ++i) {
        fx += weights[i] * abs(Dx[i]);
    }
    // kappa^T max(0, 1-c)
    if (params->opt.use_temporal_penalty) { 
        for (size_t i = 0; i < params->prob.nnodes; ++i) { 
            fx += kappa[i] * max( 0.0f, 1.0f - sol_c[i] );
        }
    }
    free(Doccx);
    free(Dx);

    return fx;
}

void solver_primaldual::print_cost() {
    //
    float cumfx = 0.0f;
    float fx = 0.0f;

    float *Doccx, *Dx;
    posix_memalign((void**)&Doccx, MEM_ALIGN_SIZE, sizeof(float)*params->prob.nocc_constraints );
    posix_memalign((void**)&Dx,    MEM_ALIGN_SIZE, sizeof(float)*params->prob.nedges );

    if (params->prob.solve_pixelwise) {
        apply_pixelwise_gradient_op(Dx, sol_c, params->prob.rows, params->prob.cols);
    } else {
        apply_sparse_op<int>(Dx, params->prob.D, sol_c );
    }
    apply_sparse_op<int>(Doccx, params->prob.Docc, sol_c);
        
    // \tau 1^T c
    fx = 0;
    for (size_t i = 0; i < params->prob.nnodes; ++i) { 
        fx += tau1[i] * sol_c[i]; 
    }
    cumfx += fx;
    printf("tau1^T c \t\t\t= %f\n", fx);

    // \lambda^T max(0, 1 - Docc x)
    fx = 0;
    for (size_t i = 0; i < params->prob.nocc_constraints; ++i) { 
        fx += occweights[i] * max( 0.0f, 1.0f - Doccx[i] );
    }
    printf("lambda^T max(0, 1-Docc c) \t= %f\n", fx);
    cumfx += fx;
    
    // |WDx|_1
    fx = 0;
    for (size_t i = 0; i < params->prob.nedges; ++i) {
        fx += weights[i] * abs(Dx[i]);
    }
    printf("|| W D c ||_1 \t\t\t= %f\n", fx);
    cumfx += fx;

    // kappa^T max(0,1-c)
    if (params->opt.use_temporal_penalty) { 
        fx = 0;
        for (size_t i = 0; i < params->prob.nnodes; ++i) {
            fx += kappa[i] * max(0.0f, 1-sol_c[i]);
        }
        printf("kappa^T max(0,1-c) \t\t= %f\n", fx);
        cumfx += fx;
    }
    printf("total \t\t\t\t= %f\n", cumfx);
    free(Doccx);
    free(Dx);    
}



void solver_primaldual::solve_y1(float* DX, float sigma_y) { 
#if defined(HAVE_SSE) || defined(HAVE_AVX)
    solve_y1_sse(DX, sigma_y);
    return;
#else
    for (size_t i = 0; i < params->prob.nedges; ++i) { 
        float y1_hat = sol_y1[i] + sigma_y*DX[i];
        int s = y1_hat>0?-1:+1;
        sol_y1[i] = s*max( -weights[i], s*y1_hat );
    }
#endif
}

#if defined(HAVE_SSE) || defined(HAVE_AVX)
void solver_primaldual::solve_y1_sse(float* DX, float sigma_y) {
    size_t n = params->prob.nedges;
    size_t nloops = n/MSTEP;

    MVAR dx, y, yhat, w, temp, abs_yhat, sign_yhat;
    MVAR m_sigma_y  = MM_set1_ps(sigma_y);
    MVAR m_zero     = MM_set1_ps(0.0f);
    MVAR m_one      = MM_set1_ps(1.0f);
    // init sign vectors:
    unsigned int ABS_MASK[MSTEP] __attribute__((aligned(32)));
    unsigned int SIGN_MASK[MSTEP] __attribute__((aligned(32)));
    for (int i = 0; i < MSTEP; ++i) { 
        ABS_MASK[i]  = 0x80000000;
        SIGN_MASK[i] = 0x7FFFFFFF;
    }
    MVARi m_signmask = MM_load_si( (MVARi*)SIGN_MASK );
    MVARi m_absmask  = MM_load_si( (MVARi*)ABS_MASK );

    for (size_t i = 0; i < nloops; ++i) {
        dx      = MM_load_ps(DX +       MSTEP*i);
        y       = MM_load_ps(sol_y1 +   MSTEP*i);
        w       = MM_load_ps(weights +  MSTEP*i);
        // yhat = y[i] + sigma*DX[i]
        yhat    = MM_add_ps(y, MM_mul_ps(m_sigma_y, dx) );  
        // sign(y) min(w, abs(y) )
        abs_yhat    = MM_andnot_ps( (MVAR)m_absmask, yhat);   
        sign_yhat   = MM_or_ps( MM_andnot_ps( (MVAR)m_signmask, yhat), m_one);   
        yhat = MM_mul_ps(sign_yhat, MM_min_ps(w, abs_yhat) );

        MM_store_ps(sol_y1 + MSTEP*i, yhat);
    }
    for (size_t i = MSTEP*nloops; i < n; ++i) {
        float y = sol_y1[i] + sigma_y*DX[i];
        float w = weights[i];
        sol_y1[i] = max(-w, min(0.0f, y)) + min(+w, max(0.0f, y));
    }
}
#endif


float solver_primaldual::solve_c(float* DTY, float* DOCCTY, float sigma, float* tau) { 
    #if defined(HAVE_SSE) || defined(HAVE_AVX)
        return solve_c_sse(DTY, DOCCTY, sigma, tau);
    #else
    // upper bound. this is a little bit faster.
    float c_prev;
    float dx_ = 0.0f;
    float dc;
    if (params->opt.layer_upper_bound > 0) { 
        for (size_t i = 0; i < params->prob.nnodes; ++i) {
            c_prev = sol_c[i];
            sol_c[i] = min(params->opt.layer_upper_bound,
                    max(0.0f, sol_c[i] - sigma*(DTY[i] + DOCCTY[i] + tau[i]) ) );
            dc = sol_c[i] - c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
        }
    } else {
        for (size_t i = 0; i < params->prob.nnodes; ++i) {
            c_prev = sol_c[i];
            sol_c[i] = max(0.0f, sol_c[i] - sigma*(DTY[i] + DOCCTY[i] + tau[i]) );
            dc = sol_c[i] - c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
        }
    }
    return dx_;
    #endif
}

float solver_primaldual::solve_c_temporal(float* DTY, 
                                    float* DOCCTY, float sigma, float* tau) { 
    #if defined(HAVE_SSE) || defined(HAVE_AVX)
        return solve_c_temporal_sse(DTY, DOCCTY, sigma, tau);
    #else
    // use a separate block for the case that hte user specified an
    // upper bound. this is a little bit faster.
    float c_prev;
    float dx_ = 0.0f;
    float dc;

    if (params->opt.layer_upper_bound > 0) {
        for (size_t i = 0; i < params->prob.nnodes; ++i) {
            float sigmatau = sigma*tau[i];
            float oneplussigmatau = 1+sigma*tau[i];
            c_prev = sol_c[i];
            // ----------------------------------------------------------------
            float z = sol_c[i] - sigma*(DTY[i]+DOCCTY[i]);
            if (z > oneplussigmatau) {
                sol_c[i] = min(z - sigmatau, params->opt.layer_upper_bound);
            } else if (z < oneplussigmatau - sigma*kappa[i]) { 
                sol_c[i] = max(0.0f, z - sigmatau + sigma*kappa[i]);
            } else {
                sol_c[i] = 1;
            }
            // ----------------------------------------------------------------
            dc = sol_c[i] - c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
        }
    } else {
        for (size_t i = 0; i < params->prob.nnodes; ++i) {
            float sigmatau = sigma*tau[i];
            float oneplussigmatau = 1+sigma*tau[i];

            c_prev = sol_c[i];
            // ----------------------------------------------------------------
            float z = sol_c[i] - sigma*(DTY[i]+DOCCTY[i]);
            if (z > oneplussigmatau) {
                sol_c[i] = z - sigmatau;
            } else if (z < oneplussigmatau - sigma*kappa[i]) { 
                sol_c[i] = max(0.0f, z - sigmatau + sigma*kappa[i]);
            } else {
                sol_c[i] = 1;
            }
            // ----------------------------------------------------------------
            dc = sol_c[i] - c_prev;
            sol_x[i] = sol_c[i] + dc;
            dx_ += dc*dc;
        }
    }
    return dx_;
    #endif
}

#if defined(HAVE_SSE) || defined(HAVE_AVX)
float solver_primaldual::solve_c_sse(float* DTY,
                            float* DOCCTY, float sigma, float* tau) {
    size_t n = params->prob.nnodes;
    size_t nloops = n/MSTEP;
    MVAR c, c_prev, x, dty, doccty, dc, dc2, temp, m_tau;
    MVAR dxsum      = MM_set1_ps(0.0f);
    MVAR m_sigma    = MM_set1_ps(sigma);
    MVAR m_zero     = MM_set1_ps(0.0f);

    if (params->opt.layer_upper_bound > 0) {
        MVAR m_ub = MM_set1_ps(params->opt.layer_upper_bound);
        for (size_t i = 0; i < nloops; i++) { 
            m_tau   = MM_load_ps(tau +      MSTEP*i);
            c       = MM_load_ps(sol_c +    MSTEP*i);
            dty     = MM_load_ps(DTY +      MSTEP*i);
            doccty  = MM_load_ps(DOCCTY +   MSTEP*i);
            c_prev  = c;
            temp    = MM_add_ps(dty, doccty);
            temp    = MM_add_ps(temp, m_tau);
            temp    = MM_mul_ps(temp, m_sigma); // sigma_c*(DTY+DOCCTY+tau)
            c       = MM_sub_ps(c, temp);
            c       = MM_max_ps(c, m_zero); // max(0, c)
            c       = MM_min_ps(c, m_ub);
            dc      = MM_sub_ps(c, c_prev);

            x       = MM_add_ps(c, dc); // x = sol_c + dc
            MM_store_ps(sol_c + MSTEP*i, c);

            MM_stream_ps(sol_x + MSTEP*i, x);
            //MM_store_ps(sol_x + MSTEP*i, x);
           
            dc2     = MM_mul_ps(dc, dc);
            dxsum   = MM_add_ps(dxsum, dc2);
        }
    } else {
        for (size_t i = 0; i < nloops; i++) { 
            m_tau   = MM_load_ps(tau +      MSTEP*i);            
            c       = MM_load_ps(sol_c +    MSTEP*i);
            dty     = MM_load_ps(DTY +      MSTEP*i);
            doccty  = MM_load_ps(DOCCTY +   MSTEP*i);
            c_prev  = c;
            temp    = MM_add_ps(dty, doccty);
            temp    = MM_add_ps(temp, m_tau);
            temp    = MM_mul_ps(temp, m_sigma); // sigma_c*(DTY+DOCCTY+tau)
            c       = MM_sub_ps(c, temp);
            c       = MM_max_ps(c, m_zero); // max(0, c)

            dc      = MM_sub_ps(c, c_prev);
            x       = MM_add_ps(c, dc); // x = sol_c + dc

            MM_store_ps(sol_c + MSTEP*i, c);
            MM_stream_ps(sol_x + MSTEP*i, x);
            //MM_store_ps(sol_x + MSTEP*i, x);
           
            dc2     = MM_mul_ps(dc, dc);
            dxsum   = MM_add_ps(dxsum, dc2);
        }
    }
    // need to verify:
    float dxsum_[MSTEP];
    MM_store_ps(dxsum_, dxsum);
    float dx = 0; for (int i = 0; i < MSTEP; ++i) { dx += dxsum_[i]; }

    for (size_t i = MSTEP*nloops; i < params->prob.nnodes; ++i) {
        float c_prev = sol_c[i];
        sol_c[i] = max(0.0f, sol_c[i] - sigma*(DTY[i] + DOCCTY[i] + tau[i]) );
        if (params->opt.layer_upper_bound > 0) {
            sol_c[i] = min(sol_c[i], params->opt.layer_upper_bound);
        }
        float dc = sol_c[i] - c_prev;
        sol_x[i] = sol_c[i] + dc;
        dx += dc*dc;
    }
    return dx;
}


float solver_primaldual::solve_c_temporal_sse(float* DTY,
                            float* DOCCTY, float sigma, float* tau) {
    size_t n = params->prob.nnodes;
    size_t nloops = n/MSTEP;
    
    float dc;

    MVAR c, m_cprev, dty, doccty, m_dc, m_dc2, temp;
    MVAR m_tau, m_sigmatau, m_oneplussigmatau;
    MVAR dxsum      = MM_set1_ps(0.0f);
    MVAR m_sigma    = MM_set1_ps(sigma);
    MVAR m_zero     = MM_set1_ps(0.0f);
    MVAR m_one      = MM_set1_ps(1.0f);
    MVAR m_ffff     = MM_cmpgt_ps( m_one, m_zero);

    MVAR m_ub       = MM_set1_ps( params->opt.layer_upper_bound );
    if (params->opt.layer_upper_bound <= 0) {
        m_ub = MM_set1_ps( 99999999.0f ); // less code.
    }
    
    MVAR m_kappa, m1, m2, m3, out1, out2, out3;

    for (size_t i = 0; i < nloops; ++i) { 
        // ----------------------------------------------------------------
        // first just compute sol_c <- sol_c - sigma*(DTY[i]+DOCCTY[i])
        m_tau               = MM_load_ps(tau + MSTEP*i);
        m_sigmatau          = MM_mul_ps(m_sigma, m_tau);
        m_oneplussigmatau   = MM_add_ps(m_one, m_sigmatau);

        c                   = MM_load_ps(sol_c +    MSTEP*i);
        m_cprev             = c;
        dty                 = MM_load_ps(DTY +      MSTEP*i);
        doccty              = MM_load_ps(DOCCTY +   MSTEP*i);
        temp                = MM_add_ps(dty, doccty);
        temp                = MM_mul_ps(temp, m_sigma); // sigma_c*(DTY+DOCCTY)
        c                   = MM_sub_ps(c, temp);
        // ----------------------------------------------------------------

        m_kappa             = MM_load_ps(kappa + MSTEP*i);
        MVAR m_sigmakappa   = MM_mul_ps(m_sigma, m_kappa);
        MVAR m_lobound      = MM_sub_ps(m_oneplussigmatau, m_sigmakappa);

        // m1 = ffff , when c >= 1+sigma*tau
        m1 = MM_cmpgt_ps(c, m_oneplussigmatau);
        // m2 = ffff , when c <= 1+sigma*tau - sigma*kappa
        m2 = MM_cmplt_ps(c, m_lobound);
        // m3 = ffff , when m1=0 & m2 = 0
        m3 = MM_xor_ps( MM_or_ps(m1, m2), m_ffff );

        // now compute all outputs:
        out1 = MM_sub_ps( c, m_sigmatau);
        out2 = MM_add_ps( out1, m_sigmakappa);
        out1 = MM_min_ps( out1, m_ub );
        out2 = MM_max_ps( out2, m_zero );
        out3 = m_one;
        // and all three, then or them:
        out1 = MM_and_ps(out1, m1);
        out2 = MM_and_ps(out2, m2);
        out3 = MM_and_ps(out3, m3);
        c = MM_or_ps( MM_or_ps(out1, out2), out3);

        MM_store_ps( sol_c + MSTEP*i, c);
        //---------------------------------------------------------------------
        m_dc = MM_sub_ps(c, m_cprev);
        temp = MM_add_ps(c, m_dc );
        m_dc2 = MM_mul_ps(m_dc, m_dc);
        dxsum = MM_add_ps(m_dc2, dxsum);
        MM_stream_ps( sol_x + MSTEP*i, temp);       
        //MM_store_ps( sol_x + MSTEP*i, temp);
    }
    float dxsum_[MSTEP];
    MM_store_ps(dxsum_, dxsum);
    float dx = 0; for (int i = 0; i < MSTEP; ++i) { dx += dxsum_[i]; }
     
    for (size_t i = MSTEP*nloops; i < params->prob.nnodes; ++i) { 
        float sol_cprev = sol_c[i];
        float sigmatau = sigma*tau[i];
        float oneplussigmatau = 1 + sigmatau;
        sol_c[i] = sol_c[i] - sigma*(DTY[i] + DOCCTY[i] );

        if (sol_c[i] > oneplussigmatau) {
            sol_c[i] = min(sol_c[i] - sigmatau, params->opt.layer_upper_bound);
        } else if (sol_c[i] < oneplussigmatau - sigma*kappa[i]) { 
            sol_c[i] = max(0.0f, sol_c[i] - sigmatau + sigma*kappa[i]);
        } else {
            sol_c[i] = 1;
        }
        dc = sol_c[i] - sol_cprev;
        sol_x[i] = sol_c[i] + dc;
        dx += dc*dc;
    }
    return dx;
}
#endif

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
    err = (params->prob.tau1.size()==0);
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
