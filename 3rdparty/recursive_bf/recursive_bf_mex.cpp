#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "qx_basic.h"
#include "qx_ppm.h"
#include "qx_recursive_bilateral_filter.h"

using namespace std;
#define throw_error(x) mexPrintf( (x) ); return;

// out = recursive_bf_mex(I, sigma_spatial=0.03, sigma_range=0.1, filtertype={0,1}, nriteration)
void print_usage() {
    mexPrintf("out = recursive_bf_mex(I, sigma_spatial=0.03, sigma_range=0.1, filtertype={0,1}, nriteration=10)\n");
    mexPrintf("\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if ( (nrhs != 5) || (nlhs != 1)) {
        mexPrintf("incorrect number of arguments.\n\n");
        print_usage();
        return;
    }
    if ((rand()%100)==1){mexPrintf("meow (o.o) meow\n");}
    float sigma_spatial = (float)mxGetScalar( prhs[1] );
    float sigma_range = (float)mxGetScalar( prhs[2] );
    int filter_type = (int)mxGetScalar( prhs[3] );
    int nr_iteration = (int)mxGetScalar( prhs[4] );
    
    int dim = mxGetNumberOfDimensions(prhs[0]);
    if (dim==2) { throw_error("image should be RGB or multiple channels for flow.\n"); }
    const mwSize* sz = mxGetDimensions(prhs[0]);
    int rows = sz[0];
    int cols = sz[1];
    int chan = sz[2];
    if (mxIsUint8(prhs[0])) {
        if (chan != 3) { throw_error("the uchar data is not an RGB image (doesnt have 3 chans).\n"); }
        
        unsigned char* img = (unsigned char*)mxGetPr( prhs[0] );
        unsigned char*** texture = qx_allocu_3(rows, cols, chan);
        double*** image = qx_allocd_3(rows, cols, chan);
        double*** img_filtered = qx_allocd_3(rows, cols, chan);
        double*** temp = qx_allocd_3(rows, cols, chan);
        double*** temp2 = qx_allocd_3(2, cols, chan);

        for(int y=0; y<rows; y++) {
            for(int x=0; x < cols; x++) {
                for(int c=0; c<chan; c++) {
                    image[y][x][c] = img[y + x*rows + c*rows*cols];
                    texture[y][x][c] = img[y + x*rows + c*rows*cols];
                }
            }
        }

        if(filter_type==0) {
            for(int i=0; i < nr_iteration; i++) {
                qx_gradient_domain_recursive_bilateral_filter(
                    img_filtered, image, texture, sigma_spatial, sigma_range, rows, cols, temp, temp2);
            }
        } else {
            double**temp_factor=qx_allocd(rows*2+2,cols);
            for(int i = 0; i < nr_iteration; i++) { 
                qx_recursive_bilateral_filter(
                    img_filtered, image, texture, sigma_spatial, sigma_range, rows, cols, temp, temp2,
                        temp_factor,&(temp_factor[rows]),&(temp_factor[rows+rows]));
            }
            qx_freed(temp_factor); temp_factor=NULL;
        }

        plhs[0] = mxCreateNumericArray(3, sz, mxUINT8_CLASS, mxREAL);
        unsigned char* out = (unsigned char*)mxGetPr(plhs[0]);
        for (int r = 0; r < rows; ++r) { 
            for (int c = 0; c < cols; ++c) {
                for (int ch = 0; ch < chan; ++ch) { 
                    out[r + c*rows + rows*cols*ch] = (unsigned char)img_filtered[r][c][ch];
                }
            }
        }
    
        qx_freeu_3(texture);
        qx_freed_3(image);
        qx_freed_3(img_filtered);
        qx_freed_3(temp);
        qx_freed_3(temp2);
    } else if (mxIsDouble(prhs[0])) {
        if (chan != 2) { throw_error("the double data doesnt look like flow (doesnt have 2 chans).\n"); }

        double* img = (double*)mxGetPr( prhs[0] );
        double*** texture = qx_allocd_3(rows, cols, 3);
        double*** image = qx_allocd_3(rows, cols, 3);
        double*** img_filtered = qx_allocd_3(rows, cols, 3);
        double*** temp = qx_allocd_3(rows, cols, 3);
        double*** temp2 = qx_allocd_3(2, cols, 3);

        for(int y=0; y<rows; y++) {
            for(int x=0; x < cols; x++) {
                for(int c=0; c<2; c++) {
                    image[y][x][c] = img[y + x*rows + c*rows*cols];
                    texture[y][x][c] = img[y + x*rows + c*rows*cols];
                }
                image[y][x][2] = img[y + x*rows + 1*rows*cols];
                texture[y][x][2] = img[y + x*rows + 1*rows*cols];
            }
        }

        double**temp_factor = qx_allocd(rows*2+2,cols);
        for (int i = 0; i < nr_iteration; i++) { 
            qx_recursive_bilateral_filter_2chan(
                   img_filtered, image, texture, sigma_spatial, sigma_range, rows, cols, temp, temp2,
                        temp_factor,&(temp_factor[rows]),&(temp_factor[rows+rows]));
        }

        qx_freed(temp_factor); temp_factor=NULL;

        plhs[0] = mxCreateNumericArray(3, sz, mxDOUBLE_CLASS, mxREAL);
        double* out = (double*)mxGetPr(plhs[0]);
        for (int r = 0; r < rows; ++r) { 
            for (int c = 0; c < cols; ++c) {
                for (int ch = 0; ch < chan; ++ch) { 
                    out[r + c*rows + rows*cols*ch] = (double)img_filtered[r][c][ch];
                }
            }
        }
        qx_freed_3(texture);
        qx_freed_3(image);
        qx_freed_3(img_filtered);
        qx_freed_3(temp);
        qx_freed_3(temp2);
    }
}