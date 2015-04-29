#include "sparse_op.h"

sparse_op<int> create_diff_op(int rows, int cols, int chan, char type) {
    sparse_op<int> D;
    int n = rows*cols*chan;
    D.R = vector<int>();
    D.C = vector<int>();
    D.val = vector<int>();
    D.R.reserve(2*n);
    D.C.reserve(2*n);
    D.val.reserve(2*n);
    D.nrows = n;
    D.ncols = n;   
 
    if ((type == 'V') || (type == 'v')) {
        for (int r = 0; r < rows; ++r) { 
            int rp = r+1;
            if (rp == rows) { rp = 0; }
            for (int c = 0; c < cols; ++c) {
                for (int k = 0; k < chan; ++k) { 
                    int i1 = linear_index(r, c, k, cols, chan);
                    int i2 = linear_index(rp, c, k, cols, chan);
                    D.R.push_back(i1);
                    D.C.push_back(i1);
                    D.val.push_back(-1);
                     
                    D.R.push_back(i1);
                    D.C.push_back(i2);
                    D.val.push_back(+1);
                }
            }
        } 
    } else if ((type == 'H') || (type == 'h')) { 
        for (int r = 0; r < rows; ++r) { 
            for (int c = 0; c < cols; ++c) { 
                int cp = c+1;
                if (cp == cols) { cp = 0; }
                for (int k = 0; k < chan; ++k) { 
                    int i1 = linear_index(r, c, k, cols, chan);
                    int i2 = linear_index(r, cp, k, cols, chan);    
                    D.R.push_back(i1);
                    D.C.push_back(i1);
                    D.val.push_back(-1);
                     
                    D.R.push_back(i1);
                    D.C.push_back(i2);
                    D.val.push_back(+1);
                }
            }
        }
    } else { 
        printf("Unrecognized diff-op type\n");
    }
    return D;
}

void create_pixelwise_diff_op(diff_op& Dx, diff_op& Dy, int rows, int cols, int chan) {
    int n = rows*cols*chan;
    // initialize
    Dx.Cpos = vector<int>();
    Dx.Cneg = vector<int>();
    Dx.Ctpos = vector<int>();
    Dx.Ctneg = vector<int>();
    
    Dy.Cpos = vector<int>();
    Dy.Cneg = vector<int>();
    Dy.Ctpos = vector<int>();
    Dy.Ctneg = vector<int>();

    Dx.Cpos.reserve(n);
    Dx.Cneg.reserve(n);
    Dx.Ctpos.reserve(n);
    Dx.Ctneg.reserve(n);

    // vertical:
    for (int r = 0; r < rows; ++r) { 
        int rn = r+1;
        int rp = r-1;
        if (rp < 0) { rp = rows-1; }
        if (rn == rows) { rn = 0; }
        for (int c = 0; c < cols; ++c) {
            for (int k = 0; k < chan; ++k) { 
                int i1 = linear_index(r, c, k, cols, chan);
                int i2 = linear_index(rn, c, k, cols, chan);
                int i0 = linear_index(rp, c, k, cols, chan);
                Dy.Cneg.push_back(i1);
                Dy.Cpos.push_back(i2);
                
                Dy.Ctneg.push_back(i1);
                Dy.Ctpos.push_back(i0);
            }
        }
    } 
    // horizontal:
    for (int r = 0; r < rows; ++r) { 
        for (int c = 0; c < cols; ++c) { 
            int cn = c+1;
            int cp = c-1;
            if (cn == cols) { cn = 0; }
            if (cp < 0) { cp = cols-1; }
            for (int k = 0; k < chan; ++k) { 
                int i1 = linear_index(r, c, k, cols, chan);
                int i2 = linear_index(r, cn, k, cols, chan);    
                int i0 = linear_index(r, cp, k, cols, chan);
                Dx.Cneg.push_back(i1);
                Dx.Cpos.push_back(i2);
                
                Dx.Ctneg.push_back(i1);
                Dx.Ctpos.push_back(i0);
            }
        }
    }
}


sparse_op<int> create_diff_op(int rows, int cols, int chan) {
    sparse_op<int> dx = create_diff_op(rows, cols, chan, 'h');
    sparse_op<int> dy = create_diff_op(rows, cols, chan, 'v');
    int n = dx.R.size(); // number on nnz elements

    sparse_op<int> d;
    d.R = vector<int>(2*n);
    d.C = vector<int>(2*n);
    d.val = vector<int>(2*n);
    d.nrows = dx.nrows*2;
    d.ncols = dx.ncols;
    // copy data:
    for (int i = 0; i < n; ++i) { 
        d.R[i] = dx.R[i];
        d.C[i] = dx.C[i];
        d.val[i] = dx.val[i];
        d.R[n+i] = dy.R[i] + dx.nrows;
        d.C[n+i] = dy.C[i];
        d.val[n+i] = dy.val[i];
    }
    return d;
}

sparse_op<int> create_diff_op_from_data(const vector<int>& p1, 
        const vector<int>& p2, int n_opcols) {
    sparse_op<int> D;
    int m = p1.size();
    D.R = vector<int>();
    D.C = vector<int>();
    D.val = vector<int>();
    D.R.reserve(2*m);
    D.C.reserve(2*m);
    D.val.reserve(2*m);
    D.nrows = m;
    D.ncols = n_opcols;   

    // each row is of the form: (Dx)_i = x[ p2[i] ] - x[ p1[i] ]
    for (int i = 0; i < m; ++i) {  
        D.R.push_back( i );
        D.C.push_back( p1[i] );
        D.val.push_back( -1 );

        D.R.push_back( i );
        D.C.push_back( p2[i] );
        D.val.push_back( +1 );
    }
    return D;
}


void apply_pixelwise_diff_op(vector<float>& out, 
        const diff_op& A, const vector<float>& in, char xpos) {
    bool transposed = ((xpos == 't') || (xpos=='T'));
    size_t m;
    if (transposed==0) { m = A.Cpos.size();
    } else { m = A.Ctpos.size(); } 
    if (out.size() != m) { out = vector<float>(m); }

    if (transposed==0) {
        for (size_t i = 0; i < m; ++i) { 
            out[i] = in[ A.Cpos[i] ] - in[ A.Cneg[i] ];
        }
    } else {
        for (size_t i = 0; i < m; ++i) { 
            out[i] = in[ A.Ctpos[i] ] - in[ A.Ctneg[i] ];
        }
    }
}


void apply_pixelwise_gradient_op(vector<float>& out, 
        const diff_op& Ax, const diff_op& Ay, const vector<float>& in) {
    size_t m;
    m = Ax.Cpos.size(); 
    if (out.size() != 2*m) { out = vector<float>(2*m); }
    
    for (size_t i = 0; i < m; ++i) { 
        out[i] = in[ Ax.Cpos[i] ] - in[ Ax.Cneg[i] ];
        out[i+m] = in[ Ay.Cpos[i] ] - in[ Ay.Cneg[i] ];
    }
}

void apply_pixelwise_gradient_op_transpose(vector<float>& out, 
        const diff_op& Ax, const diff_op& Ay, const vector<float>& in) {
    size_t m;
    m = Ax.Cpos.size(); 
    if (out.size() != m) { out = vector<float>(m); }

    for (size_t i = 0; i < m; ++i) { 
        out[i] = (in[ Ax.Ctpos[i] ]-in[ Ax.Ctneg[i] ]) + 
                 (in[ m + Ay.Ctpos[i] ] - in[m + Ay.Ctneg[i] ]);
    }
}



