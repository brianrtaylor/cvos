% mex_setup
cmd = sprintf('mex recursive_bf_mex.cpp -I. qx_recursive_bilateral_filter.cpp qx_ppm.cpp qx_basic.cpp');
eval(cmd);
