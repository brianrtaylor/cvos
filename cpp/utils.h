//! @file utils.h
//! Small general utilities.
#ifndef UTILS_H_
#define UTILS_H_

#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

typedef unsigned char uchar;
typedef unsigned long ulong;

#define ISNAN(x) ((x)!=(x))

//! time logging utility
unsigned long now_us();
double now_ms();

//! image element access operation                                              
int inline linear_index(int r, int c, int k, int cols, int chan) {              
    return k + (c + r*cols)*chan;                                               
}
//! image element access for single-channel images                                              
int inline linear_index(int r, int c, int cols) {                               
    return c + r*cols;                                                          
}                                  
//! gets image row of single channel images                                        
int inline getrow(int i, int cols) {                                            
    return i / cols;                                                            
}                                                                               
//! gets image row of single channel images                                        
int inline getcol(int i, int cols) {                                            
    return i - (i / cols)*cols;                                                 
}

//! Checks a vector for NaNs and replaces them with a given value.
template <class T> int replace_nans(vector<T>& data, T replace) {
    int n = 0;
    for (size_t i = 0; i < data.size(); ++i) { 
        if ISNAN(data[i]) { data[i] = replace; n++; }
    }
    return n;
}

//! Returns a 32x3 vector of JET colormap values
vector<uchar> init_jet_colormap();

//! Returns a 32x3 vector of HOT colormap values
vector<uchar> init_hot_colormap();

//! colormap enum for apply_colormap
enum {CMAP_HOT, CMAP_JET};

//! returns a 3-channel RGB image 
vector<uchar> apply_colormap(const vector<float>& data, int map);

//! Returns a combination of two images, weighted by alpha
template <class T> vector<uchar> overlay_image(const vector<T>& img, const vector<uchar>& map, float alpha) {
    int n = img.size();
    int m = map.size();
    vector<uchar> out = vector<uchar>(m);
    if (n==m) { // 3 channel image passed:
        for (int i = 0; i < m/3; ++i) { 
            T val = (img[3*i+0]+img[3*i+1]+img[3*i+2])/3.0;
            out[3*i+0] = alpha*val + (1-alpha)*map[3*i+0];
            out[3*i+1] = alpha*val + (1-alpha)*map[3*i+1];
            out[3*i+2] = alpha*val + (1-alpha)*map[3*i+2];
        }
    } else {
        for (int i = 0; i < m/3; ++i) { 
            out[3*i+0] = alpha*img[i] + (1-alpha)*map[3*i+0];
            out[3*i+1] = alpha*img[i] + (1-alpha)*map[3*i+1];
            out[3*i+2] = alpha*img[i] + (1-alpha)*map[3*i+2];
        }
    }
    return out;
}
#endif // UTILS_H_
