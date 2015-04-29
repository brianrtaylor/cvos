#include <math.h>
#include <matrix.h>
#include <mex.h>

#include <cstdio>
#include <cstdlib>
#include <image.h>
#include <misc.h>
#include <pnmfile.h>
#include "segment-image.h"

using namespace std;
#define throw_error(x) mexPrintf( (x) ); return;

// seg = pff_segment_mex(I, sigma=0.5, K=500, min_size=20)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
     if ( (nrhs != 4) || (nlhs != 1)) {
         mexPrintf("incorrect number of arguments.\n\n");
         return;
     }
    unsigned char* img = (unsigned char*)mxGetPr( prhs[0] );
    float sigma = (float)mxGetScalar( prhs[1] );
    float K = (float)mxGetScalar( prhs[2] );
    int min_size = (int)mxGetScalar( prhs[3] );
    
    int dim = mxGetNumberOfDimensions(prhs[0]);
    if (dim==2) { throw_error("image should be RGB.\n"); }
    const mwSize* sz = mxGetDimensions(prhs[0]);
    int height = sz[0];
    int width = sz[1];
    int chan = sz[2];
    if (chan !=3) { throw_error("image should be RGB.\n"); }
    if (!mxIsUint8(prhs[0])) { throw_error("image should be RGB char\n"); }
    
    image<rgb> *im = new image<rgb>(width, height);
    for (int r = 0; r < height; ++r) { 
        for (int c = 0; c < width; ++c) { 
            im->access[r][c].r = img[r + c*height + height*width*0];
            im->access[r][c].g = img[r + c*height + height*width*1];
            im->access[r][c].b = img[r + c*height + height*width*2];
        }
    }
    // segmentation call.
    int num_ccs; 
    image<rgb> *seg = segment_image(im, sigma, K, min_size, &num_ccs); 
    
    // --------------------------------------------------------------------
    plhs[0] = mxCreateNumericArray(3, sz, mxUINT8_CLASS, mxREAL);
    unsigned char* out = (unsigned char*)mxGetPr(plhs[0]);
    
    for (int r = 0; r < height; ++r) { 
        for (int c = 0; c < width; ++c) { 
            out[r + c*height + height*width*0] = seg->access[r][c].r;
            out[r + c*height + height*width*1] = seg->access[r][c].g;
            out[r + c*height + height*width*2] = seg->access[r][c].b;
        }
    }
    delete seg;
    delete im;
}