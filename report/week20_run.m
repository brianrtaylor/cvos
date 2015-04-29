%-------------------------------------------------------------------------
% runs some code to generate the output I want to show
%
% @param: model: appropriate model choices follow
% * fxf  : frame by frame or basic (same as ayvaci)
% * alp  : ayvaci's PAMI12 method to temporally integrate cues
% * full : our whole system with boxes and all that jazz
%
% week17_run: to run a bunch from the shell
%-------------------------------------------------------------------------
function week17_run(seqs, model, dataset)
addpath('~/projects/tao'); setup;
if ~exist('model', 'var'); model = 'full'; end;
outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.5-full'; % full
% outpath = '/pad_local/btay/projects/tao/week17/week20-cdov-r0.6-full/';
% outpath = '/pad_local/btay/projects/tao/week17/full_20141022';
%-------------------------------------------------------------------------
% parameters
%-------------------------------------------------------------------------
PKG = struct;
PKG.model = model;
PKG.LOADFLOW = true;
PKG.m = [];
PKG.relative_weights = [0.25, 0.5, 0.5];
PKG.relative_constraint_weights = [0.5, 0.5, 0.0, 0.5];
PKG.OCCPROB = 0.25;
PKG.OCCPROBLAYERTHRESH = 0.10;
PKG.MINCONSTRAINTDIST = -1.0;
PKG.CONS_PERTURB.MAXCONSTRAINTDIST = 10.0;
PKG.versiontype = sprintf('cdov_ay_71_%s', model);
PKG.CHECKPOINT = 10;
PKG.DO_CROSSBILATERALFILTERFLOW = false; % MAGIC CHANGE FOR NEW EXPERIMENTS
PKG.DO_CONSTRAINT_DIVWEIGHT = true;

PKG.VIS = 999; % prints nothing

% weights / parameters:
PKG.TAU1 = 5e-4; % 0.001;
PKG.TAU2 = 0.5;
PKG.LAMBDA = 0.5;
PKG.PAIR = 1.0;

PKG.PROB_FG = 0.0005; %0.0005
% PKG.PROB_BG = 0.0005; %0.0005 % doesn't exist yet
% PKG.PROB_FG = 0.001; %0.0005
% PKG.PROB_BG = 0.000; %0.0005

PKG.BOXHELP = true;
PKG.PROB_BOX_FG = 0.0005; % 0.005;
PKG.BOXTVWEIGHT = 0.5;
PKG.DO_BOX_WARP_WEIGHT = true;

PKG.UNARY_LAMBDA = 0.00;
PKG.DO_CONS_PERTURB = true;
PKG.DO_CONS_NOW_PERTURB = true;

% PKG.LAYER_TEMPORAL_CONSISTENCY_WEIGHT;

PKG.DO_CONSTRAINT_WARP_WEIGHT = true;
PKG.WARP_SAFESPEEDSQUARED = 50.0;
PKG.CONSTRAINT_MEMORY_FACTOR = 0.20;

% for unity prior
PKG.OCC_INSIDE_UNITY_MEMORY_FACTOR = 0.50;
PKG.UNITY_MEMORY_FACTOR = 0.20;
PKG.TV_MEMORY_FACTOR = 0.20;
PKG.FG_MEMORY_FACTOR = 0.20;

PKG.DO_UNITY_WARP_WEIGHT = true;
PKG.UNITYHELP = true;
PKG.PROB_UNITY = 0.5; % dunno where this is used lol

PKG.POST_PERTURB = false;
PKG.CAUSAL = 1;
PKG.DO_FORBACKCAUSAL = true;

if strcmp(model, 'fullbg');
  PKG.PROB_BG = 0.001;
  outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.6-fullbg';
elseif strcmp(model, 'fullnobox'); % full but no constraint perturbation 
  PKG.BOXHELP = false;
  outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.5-fullnobox';
elseif strcmp(model, 'fullnop'); % full but no constraint perturbation 
  PKG.DO_CONS_PERTURB = false; % not necessary but just in case
  PKG.DO_CONS_NOW_PERTURB = false;
  outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.5-fullnop';
elseif strcmp(model, 'fxf');
  PKG.CAUSAL = false;
  PKG.WEIGHTSHELP = false; % should be set not to help, but just in case
  PKG.UNITYHELP = false; % should be set not to help, but just in case
  PKG.BOXHELP = false; % should be set not to help, but just in case
  PKG.PROB_FG = 0; % should be set not to help, but just in case
  PKG.DO_CONS_PERTURB = false; % not necessary but just in case
  PKG.DO_CONS_NOW_PERTURB = false;
  PKG.DO_FORBACKCAUSAL = false;
  % fxf outpath: not changing this because this should be fine
  % nothing's changed except that we have to redo the first and last frame
  outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.5-fxf'; % fxf
elseif strcmp(model, 'alp');
  PKG.CAUSAL = true;
  PKG.WEIGHTSHELP = true;
  PKG.UNITYHELP = false;
  PKG.BOXHELP = false;
  PKG.PROB_FG = 0;
  PKG.DO_FORBACKCAUSAL = false;
  PKG.DO_CONS_PERTURB = false; % not necessary but just in case
  PKG.DO_CONS_NOW_PERTURB = false;
  outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.5-alp'; % alp
elseif strcmp(model, 'wfg');
  PKG.CAUSAL = true;
  PKG.WEIGHTSHELP = false;
  PKG.UNITYHELP = false;
  PKG.BOXHELP = false;
  PKG.PROB_FG = 0.0005;
  PKG.DO_FORBACKCAUSAL = false;
  PKG.DO_CONS_PERTURB = false; % not necessary but just in case
  PKG.DO_CONS_NOW_PERTURB = false;
  outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.5-wfg'; % alp
elseif strcmp(model, 'fxfperturb');
  PKG.CAUSAL = false;
  PKG.WEIGHTSHELP = false; % should be set not to help, but just in case
  PKG.UNITYHELP = false; % should be set not to help, but just in case
  PKG.BOXHELP = false; % should be set not to help, but just in case
  PKG.PROB_FG = 0; % should be set not to help, but just in case
  PKG.DO_CONS_PERTURB = false; % not necessary but just in case
  PKG.DO_CONS_NOW_PERTURB = true;
  PKG.DO_FORBACKCAUSAL = false;
  outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.6-fxfperturb'; % fxfperturb
end
% outpath = sprintf('%s-%s/', outpath, dataset);
PKG.outpath = outpath;
if ~exist(outpath, 'dir'); createRequiredFolders(outpath); end;

%----------------------------------------------------------------------------
% run
%----------------------------------------------------------------------------
for sid = seqs;
  s = num2str(sid{1});
  % try
    fprintf('========= running seq: %s ==============\n', s);
    unix(sprintf('touch ./%s.running', s));
    cdov(sid{1}, PKG);
    unix(sprintf('rm ./%s.running', s));
    fprintf('========= finished seq: %s =============\n', s);
  % catch e
  %   fprintf('========= failed seq: %s ===============\n', s);
  %   e
  % end
end
exit
end
