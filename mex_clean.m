% mex_clean
fprintf('--------------------------------------------------------------\n')
fprintf(' Accurate Fast Marching \n')
fprintf('--------------------------------------------------------------\n')
cd('3rdparty/FastMarching_version3b');
unix('rm {functions/,shortestpath/}*.mexa64 -v');
cd('../../')

fprintf('--------------------------------------------------------------\n')
fprintf(' Recursive Bilateral Filtering \n')
fprintf('--------------------------------------------------------------\n')
cd('3rdparty/recursive_bf')
unix('rm *.mexa64 -v');
cd('../../')

fprintf('--------------------------------------------------------------\n')
fprintf(' Graph Based Image Segmentation \n')
fprintf('--------------------------------------------------------------\n')
cd('3rdparty/segment')
unix('rm *.mexa64 -v');
cd('../../')

fprintf('--------------------------------------------------------------\n')
fprintf(' Solver and other utilities \n')
fprintf('--------------------------------------------------------------\n')
cd('mex/')
unix('rm *.mexa64 -v');
cd('../')

fprintf('--------------------------------------------------------------\n')
fprintf('  Postfiltering utilities\n' );
fprintf('--------------------------------------------------------------\n')
cd('cvos_aid/postfilter/')
unix('rm *.mexa64 -v');
cd('../../')
