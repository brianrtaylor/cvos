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

//! Representation of the difference operator: 
//! out( r(i) ) = in(Cpos(i)) - in(Cneg(i)) 
//! This is a special case of the sparse operator and its application should be
//! faster.
struct diff_op {
    vector<int> Cpos;
    vector<int> Cneg;
    vector<int> Ctpos;
    vector<int> Ctneg;
};

//! Creates a sparse difference matrix (pixelwise finite difference)
//! Implements a forward difference: \f$ y[i] = x[i+1] - x[i] \f$
//! with circular boundary conditions, as. \f$ y[n-1] = x[0]-x[n-1] \f$
//! @param[in] type - specifies orientation ('h' - horizontal, 'v'- vertical)
sparse_op<int> create_diff_op(int rows, int cols, int chan, char type);

//! Creates a pair of difference operators and stacks them together, as
//! \f$ D = (D_x^T, D_y^T)^T \f$
sparse_op<int> create_diff_op(int rows, int cols, int chan);

//! Creates a pair of difference operators
void create_pixelwise_diff_op(diff_op& Dx, diff_op& Dy, int rows, int cols, int chan);


//! Creates a difference matrix from provided pairs.
//! @param[in] p1, p2 - indices of pixels in (rows,cols,chan) image. Each
//! pair corresponds to a row in the created operator
sparse_op<int> create_diff_op_from_data(const vector<int>& p1, 
        const vector<int>& p2, int n_opcols);

//void create_diff_op_from_data(diff_op& D, const vector<int>& p1, 
//        const vector<int>& p2, int rows, int cols, int chan);
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
    for (size_t u = 0; u < out.size(); ++u) { out[u] = 0.0f; }
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

//! Applies a matrix Dx to a vector. The difference operator is of the form:
//!     \f$ y(i) = x(i+1) - x(i) \f$
//! and its transpose is:
//!     \f  y(i) = x(i-1) - x(i) \f$
void apply_pixelwise_diff_op(vector<float>& out, 
        const diff_op& A, const vector<float>& in, char xpos='n');

//! Applies a concatenation of Dx and Dy to a vector; formally:
//!    \f$ D  = (D_x^T, D_y^T)^T \f$
//! If input has length n, this results in output being length 2n. 
void apply_pixelwise_gradient_op(vector<float>& out, 
        const diff_op& Ax, const diff_op& Ay, const vector<float>& in);

//! Applies a transpose of D = [Dx; Dy], which results in:
//! y = D^T x = D_x^T x(1:n/2) + D_y^T x(n/2+1:end)
void apply_pixelwise_gradient_op_transpose(vector<float>& out, 
        const diff_op& Ax, const diff_op& Ay, const vector<float>& in);

#endif // SPARSE_OP_H_
