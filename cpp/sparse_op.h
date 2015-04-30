//! @file sparse_op.h
//! Representation of 'simple' (sparse) linear operators that
//! are used by the solver. The focus is on difference operators, although
//! the structure should support any sparse linear operator.
#ifndef SPARSE_OP_H_
#define SPARSE_OP_H_

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "utils.h"

using namespace std;

//! Representation of the sparse matrix A where A( R(i), C(i) )= val(i)
template <typename T> struct sparse_op {
    vector<int> R;
    vector<int> C;
    vector<T> val;
    int nrows;
    int ncols;
};

//! Creates a difference matrix from provided pairs.
//! @param[in] p1, p2 - indices of pixels in (rows,cols,chan) image. Each
//! pair corresponds to a row in the created operator
sparse_op<int> create_diff_op_from_data(const vector<int>& p1, 
        const vector<int>& p2, int n_opcols);

// ----------------------------------------------------------------------------
//                          template functions
// ----------------------------------------------------------------------------

//! Applies sparse operator to input vector.
template <class T> void apply_sparse_op( 
        vector<float>& out, const sparse_op<T>& A, const vector<float>& in, char xpos='n') {
    bool transposed = ((xpos == 't') || (xpos=='T'));
     
    if (transposed==0) {
        out = vector<float>(A.nrows, 0.0f);
    } else {
        out = vector<float>(A.ncols, 0.0f);
    } 
    int m = A.val.size();
    if (transposed==0) {
        for (int i = 0; i < m; ++i) { 
            int r = A.R[i];
            int c = A.C[i];
            T val = A.val[i];
            out[r] += float(val)*in[c];
        }
    } else {
        for (int i = 0; i < m; ++i) { 
            int r = A.R[i];
            int c = A.C[i];
            T val = A.val[i];
            out[c] += float(val)*in[r];
        }
    }
}


//! Applies sparse operator to input vector.
template <class T> void apply_sparse_op( 
        float* out, const sparse_op<T>& A, const float* in, char xpos='n') {
    bool transposed = ((xpos == 't') || (xpos=='T'));
     
    int m = A.val.size();
    if (transposed==0) {
        memset(out, 0, sizeof(float)*A.nrows );        
        for (int i = 0; i < m; ++i) { 
            int r = A.R[i];
            int c = A.C[i];
            T val = A.val[i];
            out[r] += float(val)*in[c];
        }
    } else {
        memset(out, 0, sizeof(float)*A.ncols );
        for (int i = 0; i < m; ++i) { 
            int r = A.R[i];
            int c = A.C[i];
            T val = A.val[i];
            out[c] += float(val)*in[r];
        }
    }
}

#endif // SPARSE_OP_H_
