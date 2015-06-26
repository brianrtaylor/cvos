%
BCV_PATH = '.';
% SSE_FLAGS = '-msse -msse2 -msse3 -DHAVE_SSE';
SSE_FLAGS = '';
% mex_setup
COMMON_OBJS = sprintf('%s/bcv_diff_ops.cpp', BCV_PATH);

cmd = sprintf('mex bcv_tvsegment_mex.cpp %s/tvsegment.cpp %s -I%s CXXFLAGS="-O -fPIC %s -DHAVE_MATLAB"', ...
                BCV_PATH, COMMON_OBJS, BCV_PATH, SSE_FLAGS);
eval(cmd);
