% -------------------------------------------------------------------------
% you must run this script prior to using the package
% -------------------------------------------------------------------------

CPP_SRC = '../cpp/';
VLFEAT_DIR = '/home/vasiliy/research/toolbox/VLfeat/';
% VLFEAT_DIR = '/home/btay/Source/vlfeat-0.9.18/';

VLFEAT_SRC = [VLFEAT_DIR, 'vl/'];
VLFEAT_LIB = [VLFEAT_DIR, 'bin/glnxa64/'];

% % -------------------------------------------------------------------------
% %                   delete all previous files
% % -------------------------------------------------------------------------
unix('rm *.mexa64')

% % -------------------------------------------------------------------------
% %                   build primaldual solver
% % -------------------------------------------------------------------------
% CXXFLAGS = '-DNDEBUG -O3 -fPIC -DHAVE_MATLAB -DHAVE_SSE -DHAVE_AVX -march=native -mtune=native';
% uncomment this to TURNF OFF AVX:
CXXFLAGS = '-DNDEBUG -O3 -fPIC -DHAVE_MATLAB -DHAVE_SSE -march=native -mtune=native';
% uncomment this to TURN OFF SSE+AVX:
% CXXFLAGS = '-O -fPIC -DHAVE_MATLAB';
cmd = sprintf('mex pd_wrapper.cpp -output pd_wrapper CXXFLAGS="%s" -I%s ', CXXFLAGS, CPP_SRC);
cmd = sprintf('%s %s/solver_primaldual.cpp', cmd, CPP_SRC);
cmd = sprintf('%s %s/bcv_diff_ops.cpp', cmd, CPP_SRC);
cmd = sprintf('%s %s/sparse_op.cpp', cmd, CPP_SRC);
cmd = sprintf('%s %s/utils.cpp', cmd, CPP_SRC);
eval(cmd);
fprintf('Built solver.\n');
% % -------------------------------------------------------------------------
% %                   build GMM learning (using VLFEAT)
% % -------------------------------------------------------------------------
cmd = sprintf('mex learn_constraint_gmm_mex.cpp CFLAGS="-std=c99 -O -fPIC -DVL_DISABLE_SSE2 -DVL_DISABLE_AVX" CXXFLAGS="-O -fPIC " -I%s -fpermissive ', VLFEAT_SRC);
cmd = sprintf('%s -L%s -lvl', cmd, VLFEAT_LIB);
cmd = sprintf('%s %s/generic.c ',   cmd, VLFEAT_SRC);
cmd = sprintf('%s %s/gmm.c ',       cmd, VLFEAT_SRC);
cmd = sprintf('%s %s/kmeans.c ',    cmd, VLFEAT_SRC);
cmd = sprintf('%s %s/random.c ',    cmd, VLFEAT_SRC);
cmd = sprintf('%s %s/mathop.c ',    cmd, VLFEAT_SRC);
cmd = sprintf('%s %s/kdtree.c ',    cmd, VLFEAT_SRC);
cmd = sprintf('%s %s/host.c ',      cmd, VLFEAT_SRC);
eval(cmd);

%-------------------------------------------------------------------------
%                      build GMM learning (using VLFEAT)
%-------------------------------------------------------------------------
cmd = sprintf('mex learn_bbox_gmm_mex.cpp CXXFLAGS="-O -msse2 -msse -msse3 -fPIC " -I%s ', VLFEAT_SRC);
cmd = sprintf('%s -L%s -lvl', cmd, VLFEAT_LIB);
eval(cmd);

%-------------------------------------------------------------------------
%                      build GMM evaluation
%-------------------------------------------------------------------------
fprintf('Compiling perturb_constraints_mex.cpp\n');
cmd = sprintf('mex perturb_constraints_mex.cpp CXXFLAGS="-O -fPIC" gmm_utils.cpp');
eval(cmd);
cmd = sprintf('mex eval_gmm_bboxes_mex.cpp CXXFLAGS="-O -fPIC" gmm_utils.cpp');
eval(cmd);

% compile other local MEX files:
cppfiles = dir('*.cpp');
for ii=1:length(cppfiles)   
    if strcmpi(cppfiles(ii).name, 'pd_wrapper.cpp'), continue; end
    if strcmpi(cppfiles(ii).name, 'learn_constraint_gmm_mex.cpp'), continue; end
    if strcmpi(cppfiles(ii).name, 'learn_bbox_gmm_mex.cpp'), continue; end
    if strcmpi(cppfiles(ii).name, 'gmm_utils.cpp'), continue; end
    if strcmpi(cppfiles(ii).name, 'perturb_constraints_mex.cpp'), continue; end
    if strcmpi(cppfiles(ii).name, 'eval_gmm_bboxes_mex.cpp'), continue; end

    fprintf('Compiling: %s\n', cppfiles(ii).name);
    cmd = sprintf('mex %s', cppfiles(ii).name);
    eval(cmd);
end
