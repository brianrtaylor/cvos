#include "utils.h"

unsigned long now_us() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    unsigned long ret = tv.tv_usec;
    ret += (tv.tv_sec * 1000*1000);
    return ret;
}
double now_ms() { return (now_us()/1000.0); }

vector<uchar> apply_colormap(const vector<float>& data, int map) {
    int n = data.size();
    vector<uchar> cmap;
    if (map == CMAP_JET) { cmap = init_jet_colormap(); }
    else { cmap = init_hot_colormap(); }
    

    float maxval = *max_element(data.begin(), data.end());
    float minval = *min_element(data.begin(), data.end());
    float range = maxval-minval + 1e-5f;
    int cmap_size = cmap.size()/3;
    float gap = 1.0/float(cmap_size-1);
    
    vector<uchar> out = vector<uchar>(n*3);
    for (int i = 0; i < n; ++i) { 
        float val = (data[i]-minval)/range;
        int idx = min( cmap_size-1, (int)floor(val*(cmap_size-1)));
        int idx_next = min(cmap_size-1, (int)ceil(val*(cmap_size-1)));
        float alpha = (val - idx/float(cmap_size-1))/gap;
        // RGB:
        // out[3*i+0] = cmap[3*idx+0];
        // out[3*i+1] = cmap[3*idx+1];
        // out[3*i+2] = cmap[3*idx+2];
        // BGR:
        out[3*i+0] = cmap[3*idx+2]*(1-alpha) + cmap[3*idx_next+2]*alpha;
        out[3*i+1] = cmap[3*idx+1]*(1-alpha) + cmap[3*idx_next+1]*alpha;
        out[3*i+2] = cmap[3*idx+0]*(1-alpha) + cmap[3*idx_next+0]*alpha;
    }
    return out;
}

vector<uchar> init_jet_colormap() {
    vector<uchar> colormap = vector<uchar>(3*32);
    colormap[0]=0; colormap[1]=0; colormap[2]=160; 
    colormap[3]=0; colormap[4]=0; colormap[5]=192; 
    colormap[6]=0; colormap[7]=0; colormap[8]=224; 
    colormap[9]=0; colormap[10]=0; colormap[11]=255; 
    colormap[12]=0; colormap[13]=32; colormap[14]=255; 
    colormap[15]=0; colormap[16]=64; colormap[17]=255; 
    colormap[18]=0; colormap[19]=96; colormap[20]=255; 
    colormap[21]=0; colormap[22]=128; colormap[23]=255; 
    colormap[24]=0; colormap[25]=160; colormap[26]=255; 
    colormap[27]=0; colormap[28]=192; colormap[29]=255; 
    colormap[30]=0; colormap[31]=224; colormap[32]=255; 
    colormap[33]=0; colormap[34]=255; colormap[35]=255; 
    colormap[36]=32; colormap[37]=255; colormap[38]=224; 
    colormap[39]=64; colormap[40]=255; colormap[41]=192; 
    colormap[42]=96; colormap[43]=255; colormap[44]=160; 
    colormap[45]=128; colormap[46]=255; colormap[47]=128; 
    colormap[48]=160; colormap[49]=255; colormap[50]=96; 
    colormap[51]=192; colormap[52]=255; colormap[53]=64; 
    colormap[54]=224; colormap[55]=255; colormap[56]=32; 
    colormap[57]=255; colormap[58]=255; colormap[59]=0; 
    colormap[60]=255; colormap[61]=224; colormap[62]=0; 
    colormap[63]=255; colormap[64]=192; colormap[65]=0; 
    colormap[66]=255; colormap[67]=160; colormap[68]=0; 
    colormap[69]=255; colormap[70]=128; colormap[71]=0; 
    colormap[72]=255; colormap[73]=96; colormap[74]=0; 
    colormap[75]=255; colormap[76]=64; colormap[77]=0; 
    colormap[78]=255; colormap[79]=32; colormap[80]=0; 
    colormap[81]=255; colormap[82]=0; colormap[83]=0; 
    colormap[84]=224; colormap[85]=0; colormap[86]=0; 
    colormap[87]=192; colormap[88]=0; colormap[89]=0; 
    colormap[90]=160; colormap[91]=0; colormap[92]=0; 
    colormap[93]=128; colormap[94]=0; colormap[95]=0;
    return colormap;
}

vector<uchar> init_hot_colormap() {
    vector<uchar> colormap = vector<uchar>(3*32);
    colormap[0]=21; colormap[1]=0; colormap[2]=0; 
    colormap[3]=43; colormap[4]=0; colormap[5]=0; 
    colormap[6]=64; colormap[7]=0; colormap[8]=0; 
    colormap[9]=85; colormap[10]=0; colormap[11]=0; 
    colormap[12]=107; colormap[13]=0; colormap[14]=0; 
    colormap[15]=128; colormap[16]=0; colormap[17]=0; 
    colormap[18]=149; colormap[19]=0; colormap[20]=0; 
    colormap[21]=171; colormap[22]=0; colormap[23]=0; 
    colormap[24]=192; colormap[25]=0; colormap[26]=0; 
    colormap[27]=213; colormap[28]=0; colormap[29]=0; 
    colormap[30]=235; colormap[31]=0; colormap[32]=0; 
    colormap[33]=255; colormap[34]=0; colormap[35]=0; 
    colormap[36]=255; colormap[37]=21; colormap[38]=0; 
    colormap[39]=255; colormap[40]=43; colormap[41]=0; 
    colormap[42]=255; colormap[43]=64; colormap[44]=0; 
    colormap[45]=255; colormap[46]=85; colormap[47]=0; 
    colormap[48]=255; colormap[49]=107; colormap[50]=0; 
    colormap[51]=255; colormap[52]=128; colormap[53]=0; 
    colormap[54]=255; colormap[55]=149; colormap[56]=0; 
    colormap[57]=255; colormap[58]=171; colormap[59]=0; 
    colormap[60]=255; colormap[61]=192; colormap[62]=0; 
    colormap[63]=255; colormap[64]=213; colormap[65]=0; 
    colormap[66]=255; colormap[67]=235; colormap[68]=0; 
    colormap[69]=255; colormap[70]=255; colormap[71]=0; 
    colormap[72]=255; colormap[73]=255; colormap[74]=32; 
    colormap[75]=255; colormap[76]=255; colormap[77]=64; 
    colormap[78]=255; colormap[79]=255; colormap[80]=96; 
    colormap[81]=255; colormap[82]=255; colormap[83]=128; 
    colormap[84]=255; colormap[85]=255; colormap[86]=160; 
    colormap[87]=255; colormap[88]=255; colormap[89]=192; 
    colormap[90]=255; colormap[91]=255; colormap[92]=224; 
    colormap[93]=255; colormap[94]=255; colormap[95]=255; 
    return colormap;
}
