//! @file bcv_diff_ops.cpp
#include "bcv_diff_ops.h"

void apply_pixelwise_gradient_op(vector<float>& out, const vector<float>& in, int rows, int cols, int chan) { 
    size_t m = in.size();
    if (out.size() != 2*m) { out = vector<float>(2*m); }
    #if defined(HAVE_MATLAB) & !defined(HAVE_SSE)
    apply_dx_matlab( &out[0], &in[0], rows, cols, chan);
    apply_dy_matlab( &out[m], &in[0], rows, cols, chan);
    #elif defined(HAVE_MATLAB) & defined(HAVE_SSE)
    apply_dx_matlab_sse( &out[0], &in[0], rows, cols, chan);
    apply_dy_matlab_sse( &out[m], &in[0], rows, cols, chan);
    #else
    // TODO: implement SSE versions..
    apply_dx( &out[0], &in[0], rows, cols, chan);
    apply_dy( &out[m], &in[0], rows, cols, chan);    
    #endif
}

void apply_pixelwise_gradient_op_transpose(vector<float>& out, const vector<float>& in, int rows, int cols, int chan) { 
    size_t m = in.size();
    if (out.size() != m/2) { out = vector<float>(m/2); }
    // order of application matters!!
    #if defined(HAVE_MATLAB) & !defined(HAVE_SSE)
    apply_dxt_matlab( &out[0], &in[0], rows, cols, chan);
    apply_dyt_matlab( &out[0], &in[m/2], rows, cols, chan);
    #elif defined(HAVE_MATLAB) & defined(HAVE_SSE)
    apply_dxt_matlab_sse( &out[0], &in[0], rows, cols, chan);
    apply_dyt_matlab_sse( &out[0], &in[m/2], rows, cols, chan);            
    #else
    // TODO: implement SSE version..    
    apply_dxt( &out[0], &in[0], rows, cols, chan);
    apply_dyt( &out[0], &in[m/2], rows, cols, chan);            
    #endif
}

void apply_pixelwise_gradient_op(vector<float>& out, const vector<float>& in, int rows, int cols) { 
    size_t m = in.size();
    if (out.size() != 2*m) { out = vector<float>(2*m); }
    #if defined(HAVE_MATLAB) & !defined(HAVE_SSE)
    apply_dx_matlab( &out[0], &in[0], rows, cols);
    apply_dy_matlab( &out[m], &in[0], rows, cols);
    #elif defined(HAVE_MATLAB) & defined(HAVE_SSE)
    apply_dx_matlab_sse( &out[0], &in[0], rows, cols);
    apply_dy_matlab_sse( &out[m], &in[0], rows, cols);
    #elif defined(HAVE_SSE)
    apply_dx_sse( &out[0], &in[0], rows, cols);
    apply_dy_sse( &out[m], &in[0], rows, cols);    
    #else
    apply_dx( &out[0], &in[0], rows, cols);
    apply_dy( &out[m], &in[0], rows, cols);        
    #endif
}

void apply_pixelwise_gradient_op_transpose(vector<float>& out, const vector<float>& in, int rows, int cols) { 
    size_t m = in.size();
    if (out.size() != m/2) { out = vector<float>(m/2); }
    // order of application matters!!
    #if defined(HAVE_MATLAB) & !defined(HAVE_SSE)
    apply_dxt_matlab( &out[0], &in[0], rows, cols);
    apply_dyt_matlab( &out[0], &in[m/2], rows, cols);
    #elif defined(HAVE_MATLAB) & defined(HAVE_SSE)
    apply_dxt_matlab_sse( &out[0], &in[0], rows, cols);
    apply_dyt_matlab_sse( &out[0], &in[m/2], rows, cols);            
    #elif defined(HAVE_SSE)   
    apply_dxt_sse( &out[0], &in[0], rows, cols);
    apply_dyt_sse( &out[0], &in[m/2], rows, cols);            
    #else
    apply_dxt( &out[0], &in[0], rows, cols);
    apply_dyt( &out[0], &in[m/2], rows, cols);            
    #endif
}

void apply_dy(float* out, const float* in, int rows, int cols, int chan) { 
    for (int r = 0; r < rows-1; ++r) { 
        int rn = r+1;
        for (int c = 0; c < cols; ++c) { 
            for (int k = 0; k < chan; ++k) { 
                int i1 = linear_index(r, c, k, cols, chan);
                int i2 = linear_index(rn, c, k, cols, chan);
                out[i1] = in[i2] - in[i1];
            }
        }
    }
    //
    int r = rows-1;
    int rn = 0;
    for (int c = 0; c < cols; ++c) {
        for (int k = 0; k < chan; ++k) { 
            int i1 = linear_index(r, c, k, cols, chan);
            int i2 = linear_index(rn, c, k, cols, chan);
            out[i1] = in[i2] - in[i1];
        }
    }
}

void apply_dy(float* out, const float* in, int rows, int cols) { 
    for (int r = 0; r < rows-1; ++r) { 
        int rn = r+1;
        for (int c = 0; c < cols; ++c) { 
            int i1 = r*cols + c;
            int i2 = rn*cols + c;
            out[i1] = in[i2] - in[i1];
        }
    }
    int r = rows-1;
    int rn = 0;
    for (int c = 0; c < cols; ++c) {
        int i1 = r*cols + c;
        int i2 = rn*cols + c;
        out[i1] = in[i2] - in[i1];
    }
}

void apply_dyt(float* out, const float* in, int rows, int cols, int chan) { 
    for (int r = 1; r < rows; ++r) { 
        int rp = r-1;
        for (int c = 0; c < cols; ++c) { 
            for (int k = 0; k < chan; ++k) { 
                int i1 = linear_index(r, c, k, cols, chan);
                int i2 = linear_index(rp, c, k, cols, chan);
                //! NOTICE THAT THIS IS '+=', rather than '='.
                //! YOU MUST RUN apply_dxt prior to running this function
                out[i1] += (in[i2] - in[i1]);
            }
        }
    }
    //
    int r = 0;
    int rp = rows-1;
    for (int c = 0; c < cols; ++c) { 
        for (int k = 0; k < chan; ++k) { 
            int i1 = linear_index(r, c, k, cols, chan);
            int i2 = linear_index(rp, c, k, cols, chan);
            //! NOTICE THAT THIS IS '+=', rather than '='.
            //! YOU MUST RUN apply_dxt prior to running this function
            out[i1] += (in[i2] - in[i1]);
        }
    }
}

void apply_dyt(float* out, const float* in, int rows, int cols) { 
    for (int r = 1; r < rows; ++r) { 
        int rp = r-1;
        for (int c = 0; c < cols; ++c) { 
            int i1 = r*cols + c;
            int i0 = rp*cols + c;
            //! NOTICE THAT THIS IS '+=', rather than '='.
            //! YOU MUST RUN apply_dxt prior to running this function
            out[i1] += (in[i0] - in[i1]);
        }
    }
    int r = 0;
    int rp = rows-1;
    for (int c = 0; c < cols; ++c) { 
        int i1 = r*cols + c;
        int i0 = rp*cols + c;
        //! NOTICE THAT THIS IS '+=', rather than '='.
        //! YOU MUST RUN apply_dxt prior to running this function
        out[i1] += (in[i0] - in[i1]);
    }
}

void apply_dx(float* out, const float* in, int rows, int cols, int chan) { 
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols-1; ++c) {
            int cn = c + 1;
            for (int k = 0; k < chan; ++k) {
                int i1 = linear_index(r, c, k, cols, chan);
                int i2 = linear_index(r, cn, k, cols, chan);
                out[i1] = in[i2] - in[i1];
            }
        }
        //
        int c = cols-1;
        int cn = 0;
        for (int k = 0; k < chan; ++k) { 
            int i1 = linear_index(r, c, k, cols, chan);
            int i2 = linear_index(r, cn, k, cols, chan);
            out[i1] = in[i2] - in[i1];
        }
    }
}

void apply_dx(float* out, const float* in, int rows, int cols) { 
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols-1; ++c) { 
            int i1 = r*cols + c;
            int i2 = i1 + 1;
            out[i1] = in[i2] - in[i1];
        }
        //
        int c = cols-1;
        int cn = 0;
        int i1 = r*cols + c;
        int i2 = r*cols + cn;
        out[i1] = in[i2] - in[i1];
    }
}

void apply_dxt(float* out, const float* in, int rows, int cols, int chan) {
    for (int r = 0; r < rows; ++r) { 
        int c = 0;
        int cp = cols-1;
        for (int k = 0; k < chan; ++k) { 
            int i1 = linear_index(r, c, k, cols, chan);
            int i2 = linear_index(r, cp, k, cols, chan);
            out[i1] = (in[i2] - in[i1]);
        }
        for (int c = 1; c < cols; ++c) {
            cp = c - 1;
            for (int k = 0; k < chan; ++k) { 
                int i1 = linear_index(r, c, k, cols, chan);
                int i2 = linear_index(r, cp, k, cols, chan);
                out[i1] = (in[i2] - in[i1]);
            }
        }
    }
}

void apply_dxt(float* out, const float* in, int rows, int cols) {
    for (int r = 0; r < rows; ++r) { 
        int c = 0;
        int cp = cols-1;
        int i1 = r*cols + c;
        int i0 = r*cols + cp;
        out[i1] = (in[i0] - in[i1]);

        for (int c = 1; c < cols; ++c) {
            int i1 = r*cols + c;
            int i0 = i1 - 1;
            out[i1] = (in[i0] - in[i1]);
        }
    }
}

// ----------------------------------------------------------------------------
#ifdef HAVE_SSE
void apply_dy_sse(float* out, const float* in, int rows, int cols) { 
    __m128 m_out, m_in1, m_in2;
    int nloops_c = cols/4;
    for (int r = 0; r < rows; ++r) { 
        int rn = r+1;
        if (rn == rows) { rn = 0; }
        for (int c = 0; c < nloops_c; ++c) {
            m_in1 = _mm_loadu_ps( in + (r*cols+4*c) );
            m_in2 = _mm_loadu_ps( in + (rn*cols+4*c) );
            m_out = _mm_sub_ps(m_in2, m_in1);
            _mm_storeu_ps( out + (r*cols+4*c), m_out );
        }
        for (int c = nloops_c*4; c < cols; ++c) { 
            int i1 = r*cols + c;
            int i2 = rn*cols + c;
            out[i1] = in[i2] - in[i1];
        }
    }
}

void apply_dyt_sse(float* out, const float* in, int rows, int cols) { 
    __m128 m_out, m_in1, m_in0, m_temp;
    int nloops_c = cols/4;
    for (int r = 0; r < rows; ++r) { 
        int rp = r-1;
        if (rp < 0) { rp = rows-1; }
        for (int c = 0; c < nloops_c; ++c) {
            m_out = _mm_loadu_ps( out + (r*cols+4*c) );
            m_in1 = _mm_loadu_ps( in + (r*cols+4*c) );
            m_in0 = _mm_loadu_ps( in + (rp*cols+4*c) );
            m_temp = _mm_sub_ps(m_in0, m_in1);
            m_out = _mm_add_ps(m_out, m_temp);
            _mm_storeu_ps( out + (r*cols+4*c), m_out );
            //! NOTICE THAT THIS IS '+=', rather than '='.
        }
        for (int c = nloops_c*4; c < cols; ++c) { 
            int i1 = r*cols + c;
            int i0 = rp*cols + c;
            //! NOTICE THAT THIS IS '+=', rather than '='.
            //! YOU MUST RUN apply_dxt prior to running this function
            out[i1] += (in[i0] - in[i1]);
        }
    }
}

void apply_dx_sse(float* out, const float* in, int rows, int cols) { 
    __m128 m_out, m_in1, m_in2;
    int nloops_c = (cols-1)/4;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < (cols-1)/4; ++c) { 
            m_in1 = _mm_loadu_ps( in + (r*cols+4*c) );
            m_in2 = _mm_loadu_ps( in + (r*cols+4*c+1) );
            m_out = _mm_sub_ps(m_in2, m_in1);
            _mm_storeu_ps( out + (r*cols+4*c), m_out );
        }
        for (int c = nloops_c*4; c < cols-1; ++c) {
            int i1 = r*cols + c;
            int i2 = i1 + 1;
            out[i1] = in[i2] - in[i1];
        }
        //
        int c = cols-1;
        int cn = 0;
        int i1 = r*cols + c;
        int i2 = r*cols + cn;
        out[i1] = in[i2] - in[i1];
    }
}

void apply_dxt_sse(float* out, const float* in, int rows, int cols) { 
    __m128 m_out, m_in1, m_in0;
    int nloops_c = (cols)/4;

    for (int r = 0; r < rows; ++r) { 
        int c = 0;
        int cp = cols-1;
        int i1 = r*cols + c;
        int i0 = r*cols + cp;
        out[i1] = (in[i0] - in[i1]);

        for (int c = 1; c < 4; ++c) {
            int i1 = r*cols + c;
            int i0 = i1 - 1;
            out[i1] = (in[i0] - in[i1]);
        }
        for (int c = 1; c < nloops_c; ++c) {
            m_in1 = _mm_loadu_ps( in + (r*cols+4*c) );
            m_in0 = _mm_loadu_ps( in + (r*cols+4*c-1) );
            m_out = _mm_sub_ps(m_in0, m_in1);
            _mm_storeu_ps( out + (r*cols+4*c), m_out );
        }
        for (int c = nloops_c*4; c < cols; ++c) {
            int i1 = r*cols + c;
            int i0 = i1 - 1;
            out[i1] = (in[i0] - in[i1]);
        }
    }
}
#endif

// ----------------------------------------------------------------------------
#if defined(HAVE_MATLAB) && !defined(HAVE_SSE)

void apply_dy_matlab(float* out, const float* in, int rows, int cols, int chan) { 
    for (int k = 0; k < chan; ++k) { 
        for (int c = 0; c < cols; ++c) { 
            for (int r = 0; r < rows; ++r) { 
                int rn = r+1;
                if (r == rows-1) { rn = 0; }
                int i1 = r + c*rows + k*rows*cols;
                int i2 = rn + c*rows + k*rows*cols;
                out[i1] = in[i2] - in[i1];
            }
        }
    }
}

void apply_dyt_matlab(float* out, const float* in, int rows, int cols, int chan) { 
    for (int k = 0; k < chan; ++k) { 
        for (int c = 0; c < cols; ++c) { 
            for (int r = 0; r < rows; ++r) { 
                int rp = r-1;
                if (r == 0) { rp = rows-1; }
                int i1 = r + c*rows + k*rows*cols;
                int i0 = rp + c*rows + k*rows*cols;
                //! NOTICE THAT THIS IS '+=', rather than '='.
                //! YOU MUST RUN apply_dxt prior to running this function
                out[i1] += (in[i0] - in[i1]);
            }
        }
    }
}

void apply_dx_matlab(float* out, const float* in, int rows, int cols, int chan) { 
    for (int k = 0; k < chan; ++k) { 
        for (int c = 0; c < cols; ++c) { 
            int cn = c+1;
            if (c == cols-1) { cn = 0; }
            for (int r = 0; r < rows; ++r) {
                int i1 = r + c*rows + k*cols*rows;
                int i2 = r + cn*rows + k*cols*rows;
                out[i1] = in[i2] - in[i1];
            }
        }
    }
}

void apply_dxt_matlab(float* out, const float* in, int rows, int cols, int chan) {
    for (int k = 0; k < chan; ++k) { 
        for (int c = 0; c < cols; ++c) {
            int cp = c-1;
            if (c == 0) { cp = cols-1; }
            for (int r = 0; r < rows; ++r) { 
                int i1 = r + c*rows + k*rows*cols;
                int i0 = r + cp*rows + k*rows*cols;
                out[i1] = (in[i0] - in[i1]);
            }
        }
    }
}
#endif

// ----------------------------------------------------------------------------

#if defined(HAVE_MATLAB) && defined(HAVE_SSE)
void apply_dy_matlab_sse(float* out, const float* in, int rows, int cols, int chan) { 
    __m128 m_out, m_in1, m_in2;
    int nloops_r = (rows-4)/4;
    for (int k = 0; k < chan; ++k) { 
        for (int c = 0; c < cols; ++c) { 
            //sse loop:
            for (int r = 0; r < nloops_r; ++r) { 
                m_in1 = _mm_loadu_ps( in + (4*r+c*rows+k*rows*cols) );
                m_in2 = _mm_loadu_ps( in + (1+4*r+c*rows+k*rows*cols) );
                m_out = _mm_sub_ps(m_in2, m_in1);
                _mm_storeu_ps( out + (4*r+c*rows+k*rows*cols), m_out );
            }
            // remainder: 
            for (int r = nloops_r*4; r < rows; ++r) { 
                int rn = r+1;
                if (r == rows-1) { rn = 0; }
                int i1 = r + c*rows+k*rows*cols;
                int i2 = rn + c*rows+k*rows*cols;
                out[i1] = in[i2] - in[i1];
            }
        }
    }
}

void apply_dyt_matlab_sse(float* out, const float* in, int rows, int cols, int chan) { 
    __m128 m_out, m_in1, m_in0, m_temp;
    int nloops_r = (rows-4)/4;
    for (int k = 0; k < chan; ++k) { 
        for (int c = 0; c < cols; ++c) { 
            // first block:
            for (int r = 0; r < 4; ++r) { 
                int rp = r-1;
                if (r == 0) { rp = rows-1; }
                int i1 = r + c*rows + k*rows*cols;
                int i0 = rp + c*rows + k*rows*cols;
                //! NOTICE THAT THIS IS '+=', rather than '='.
                //! YOU MUST RUN apply_dxt prior to running this function
                out[i1] += (in[i0] - in[i1]);
            }
            // sse loop:
            for (int r = 1; r < nloops_r; ++r) { 
                m_out = _mm_loadu_ps( out + (4*r+c*rows+k*rows*cols) );
                m_in1 = _mm_loadu_ps( in + (4*r+c*rows+k*rows*cols) );
                m_in0 = _mm_loadu_ps( in + (4*r+c*rows-1+k*rows*cols) );
                m_temp = _mm_sub_ps(m_in0, m_in1);
                m_out = _mm_add_ps(m_out, m_temp);
                _mm_storeu_ps( out + (4*r+c*rows+k*rows*cols), m_out );
                //! NOTICE THAT THIS IS '+=', rather than '='.
            }
            for (int r = 4*nloops_r; r < rows; ++r) { 
                int rp = r-1;
                int i1 = r + c*rows;
                int i0 = rp + c*rows;
                //! NOTICE THAT THIS IS '+=', rather than '='.
                //! YOU MUST RUN apply_dxt prior to running this function
                out[i1] += (in[i0] - in[i1]);
            }
        }
    }
}

void apply_dx_matlab_sse(float* out, const float* in, int rows, int cols, int chan) { 
    __m128 m_out, m_in1, m_in2;
    int nloops_r = rows/4;
    for (int k = 0; k < chan; ++k) { 
        for (int c = 0; c < cols; ++c) { 
            int cn = c+1;
            if (c == cols-1) { cn = 0; }
            // sse loop:
            for (int r = 0; r < nloops_r; ++r) { 
                m_in1 = _mm_loadu_ps( in + (4*r+c*rows+k*rows*cols) );
                m_in2 = _mm_loadu_ps( in + (4*r+cn*rows+k*rows*cols) );
                m_out = _mm_sub_ps(m_in2, m_in1);
                _mm_storeu_ps( out + (4*r+c*rows+k*rows*cols), m_out );
            }
            // remainder:
            for (int r = nloops_r*4; r < rows; ++r) {
                int i1 = r + c*rows+k*rows*cols;
                int i2 = r + cn*rows+k*rows*cols;
                out[i1] = in[i2] - in[i1];
            }
        }
    }
}

void apply_dxt_matlab_sse(float* out, const float* in, int rows, int cols, int chan) {
    __m128 m_out, m_in1, m_in0;
    int nloops_r = rows/4;
    for (int k = 0; k < chan; ++k) { 
        for (int c = 0; c < cols; ++c) {
            int cp = c-1;
            if (c == 0) { cp = cols-1; }
            // sse loop:
            for (int r = 0; r < nloops_r; ++r) { 
                m_in1 = _mm_loadu_ps( in + (4*r+c*rows+k*rows*cols) );
                m_in0 = _mm_loadu_ps( in + (4*r+cp*rows+k*rows*cols) );
                m_out = _mm_sub_ps(m_in0, m_in1);
                _mm_storeu_ps( out + (4*r+c*rows+k*rows*cols), m_out );
            }
            for (int r = nloops_r*4; r < rows; ++r) { 
                int i1 = r + c*rows + k*rows*cols;
                int i0 = r + cp*rows + k*rows*cols;
                out[i1] = (in[i0] - in[i1]);
            }
        }
    }
}

#endif
// -----------------------------------------------------------------

