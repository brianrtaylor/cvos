//! cvos_common.h
#define C2MAT(x) ((x)+1)
#define MAT2C(x) ((x)-1)

int inline linear_index(int r, int c, int k, int rows, int cols) {
  return r + c*rows + k*rows*cols;
}

int inline linear_index(int r, int c, int rows) {
  return r + c*rows;
}

int inline getrow(int i, int rows) { 
    return i - (i / rows)*rows;
}

int inline getcol(int i, int rows) { 
    return (i / rows);
}