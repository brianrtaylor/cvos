% mex_setup
% you must run this script prior to using the package
fprintf('--------------------------------------------------------------\n');
fprintf(' Piotr Computer Vision Matlab Toolbox \n');
fprintf('--------------------------------------------------------------\n');
cd('3rdparty/piotr_toolbox/external');
run('toolboxCompile');
cd('../../../');

fprintf('--------------------------------------------------------------\n');
fprintf(' Accurate Fast Marching \n');
fprintf('--------------------------------------------------------------\n');
cd('3rdparty/FastMarching_version3b');
run('compile_c_files');
cd('../../');

fprintf('--------------------------------------------------------------\n');
fprintf(' Recursive Bilateral Filtering \n');
fprintf('--------------------------------------------------------------\n');
cd('3rdparty/recursive_bf');
run('mex_setup');
cd('../../');

fprintf('--------------------------------------------------------------\n');
fprintf(' Graph Based Image Segmentation \n');
fprintf('--------------------------------------------------------------\n');
cd('3rdparty/segment');
run('mex_setup');
cd('../../');

fprintf('--------------------------------------------------------------\n');
fprintf(' Solver and other utilities \n');
fprintf('--------------------------------------------------------------\n');
cd('mex/');
run('mex_setup');
cd('../');

fprintf('--------------------------------------------------------------\n');
fprintf(' Postfiltering utilities\n');
fprintf('--------------------------------------------------------------\n');
cd('cvos_aid/postfilter/');
run('mex_setup');
cd('../../');

fprintf('Appears to have succeeded.\n');
