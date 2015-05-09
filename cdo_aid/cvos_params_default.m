function params = cvos_params_default

params = struct;

params.TINY_UV = 1.0;

params.CAUSAL = true;
params.PROB_FG = 1e-4;
params.PROB_BOX_FG = 5e-4;

params.relative_weights = [0.25 0.5 0.5];
params.relative_weights_constraints = [0.25 0.5 0.0, 0.5];
params.SIGMOID_SCALE = 100;
params.SIGMOID_MEAN = 0.55;

params.WEIGHTSHELP = true;
params.WEIGHTS_LOW_CUTOFF = 0.1;
params.CONSTRAINTS_LOW_CUTOFF = 0.01;
params.E_SIGMOID_MEAN = 0.75;
params.E_SIGMOID_SCALE = 10.0;

params.DIVFLAG = false;
params.FERRARIFLAG = true;
params.LOCALFLOWADJUST = false; % not great (try with more than cube)
params.LOCALFLOWWEIGHTADJUST = false; % unused?
params.WEIGHT_MAG_THRESHOLD = 0.65;
params.THETA = 0.2;
params.MINWUV = 5e-3;
params.MAXUVADJUSTPROPORTION = 0.9;

params.DIV_UV_THRESHOLD = 0.01;

params.PAIR = 1.0;
params.LAMBDA = 0.5;
params.TAU1 = 5e-4;
params.TAU2 = 0.5;

% essentially turned off for now
params.UNARY_LAMBDA = 0.00; % 10 x stronger than PROB_FG

params.OCCPROBLAYERTHRESH = 0.1;
params.OCCPROB = 0.5;
params.CONSTRAINT_WEIGHT_THRESH = 0.10;
params.DO_CONSTRAINT_DIVWEIGHT = true;
params.MINCONSTRAINTDIST = 2.0;
params.DO_CROSSBILATERALFILTERFLOW = true;

params.WARP_SAFESPEEDSQUARED = 50.0;
params.DO_CONSTRAINT_WARP_WEIGHT = true; 
params.DO_CONSTRAINT_OBJECT_WARP_WEIGHT = true;
params.DO_FOREGROUND_WARP_WEIGHT = true;
params.DO_WEIGHT_WARP_WEIGHT = true;
% params.TV_WARP_WEIGHT = true;

params.VIS = 500;
params.LOADFLOW = true;
params.m = [];

% TODO: increase or decrease memory (weights: seems like increase for 
%   moseg, decrease for videosegmentation)
% killed if not a part of any segmentation for 5 frames
params.UNARY_CONSTRAINT_MEMORY_FACTOR = 0.2;
% killed if not a part of any segmentation for 5 frames just about
params.CONSTRAINT_MEMORY_FACTOR = 0.20; 
% killed if not a part of any segmentation for 5 frames just about
params.WEIGHT_MEMORY_FACTOR = 0.20;
% killed if not a part of any segmentation for 10 frames just about
params.FG_MEMORY_FACTOR = 0.20;
% killed if not a part of any segmentation for 5 frames just about
params.TV_MEMORY_FACTOR = 0.20;
params.CONSTRAINT_INSIDE_FG_MEMORY_FACTOR = 0.5; % kill em off in 2 frames
params.OCC_INSIDE_FG_MEMORY_FACTOR = 0.5; % kill em off in 2 frames
params.OCC_INSIDE_WEIGHTS_MEMORY_FACTOR = 0.5; % kill em off in 2 frames
% params.TV_INSIDE_FG_MEMORY_FACTOR = 0.5; % kill em off in 2 frames 
% isn't used, the weighst additing naturally does this % isn't used, 
% the weighst additing naturally does this
params.OCCMETHOD = 'RGB_RES';

% increases weight on foreground constraint each segmentation
params.RUNNING_FG_PRIOR = true;
% increases tv penalty on regions to keep same label with each segmentation
params.RUNNING_TV_PRIOR = true;

params.dsk = strel('disk', 5);

% vasiliy gmm peturb weights
params.DO_CONS_NOW_PERTURB = true;
params.DO_CONS_PERTURB = true;
params.CONS_PERTURB = struct;
params.CONS_PERTURB.GROUP_SIZE = 4;
params.CONS_PERTURB.INVALID_COST = 50; % TODO: ask vasiliy what this should be
params.CONS_PERTURB.NUM_GMM_REPETITIONS = 1;
params.CONS_PERTURB.NUM_GMM_ITERATIONS = 5; % number of GMM iterations
params.CONS_PERTURB.NUM_GMM_CLUSTERS = 3; % number of components in GMM
params.CONS_PERTURB.MAX_SHIFT = 10; % largest translation in pixels (in one direction)
params.CONS_PERTURB.NUM_ROTATION = 3; % number of rotations (should be odd)
params.CONS_PERTURB.MAX_STRETCH = 5; % largest stretch in pixels
params.CONS_PERTURB.MAX_ROTATION = pi/6; % largest rotation
params.CONS_PERTURB.ROTATION_COST = 5*0.5; % cost (per radian)
params.CONS_PERTURB.TRANSLATION_COST = 5*0.5; % cost of translation (per pixel)
params.CONS_PERTURB.EDGE_DISTANCE_COST = 50; % cost of being close to an edge
params.CONS_PERTURB.STRETCH_COST = 5*1.0; % cost of stretching
params.CONS_PERTURB.MAXCONSTRAINTDIST = 20.0; % max distance between occd and occr
params.CONS_PERTURB.DO_CONS_PERTURB = params.DO_CONS_PERTURB;

params.DO_UNARY_CONSTRAINT_WARP_WEIGHT = true;
params.DO_CONSTRAINT_WEIGHT_SIGMOID = false;

params.DO_FORBACKCAUSAL = false;

params.DO_RECOMPUTE_OLD_WEIGHTS_IN_CURRENT_FRAME = false;

% unary occ
params.UNARYOCC = false;

% boxes for edges
params.BOX_RAD = 15;
params.BOX_CONF_MEMORY_FACTOR = 0.2;
params.BOX_CONF_INSIDE_MEMORY_FACTOR = 0.5;
params.BOXTVWEIGHT = 0.5;
params.DO_BOX_WARP_WEIGHT = true;
params.BOXHELP = true;

% weights propagation
params.WEIGHTHELP = true;

% constraints propagation
params.CONSTRAINTHELP = true;

params.DO_UNITY_WARP_WEIGHT = true;
params.UNITYHELP = true;
params.PROB_UNITY = 0.5;

params.OCC_INSIDE_UNITY_MEMORY_FACTOR = 0.50;
params.UNITY_MEMORY_FACTOR = 0.20;
params.POST_PERTURB = false;

params.CHECKPOINT = 20;

params.versiontype = 'sup';
end
