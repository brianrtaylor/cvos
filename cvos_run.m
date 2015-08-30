%-----------------------------------------------------------------------------
% cvos_run
%
% runs the cvos framework on a specific sequence specified by dataPaths. The
% models listed below include the experiments shown in the paper
%
% @param: model: appropriate model choices follow
% * fxfperturb: frame by frame or basic + perturbation
% * ayv: method from Ayvaci et al. PAMI12 to temporally integrate cues
% * full: our full system 
% * fullnopcbf: full system but no flow extrapolation
%-----------------------------------------------------------------------------
function cvos_run(seqs, model, dataset)
setup;
if ~exist('model', 'var'); model = 'full'; end;

plotbase = './';
outpath = './results/';

% default parameters
PKG = struct;
PKG.model = model;
PKG.LOADFLOW = true;
PKG.m = [];
PKG.relative_weights = [0.25, 0.5, 0.5];
PKG.relative_constraint_weights = [0.5, 0.5, 0.0, 0.5];
PKG.OCCPROB = 0.05;
PKG.OCCPROBLAYERTHRESH = 0.10;
PKG.MINCONSTRAINTDIST = -1.0;
PKG.CONS_PERTURB.MAXCONSTRAINTDIST = 10.0;
PKG.versiontype = sprintf('cvos_%s', model);
PKG.CHECKPOINT = 10;
PKG.DO_CROSSBILATERALFILTERFLOW = true;
PKG.DO_CONSTRAINT_DIVWEIGHT = true;
PKG.TAU1 = 5e-4;
PKG.LAMBDA = 0.5;
PKG.PAIR = 1.0;
PKG.PROB_FG = 0.0005;
PKG.BOXHELP = true;
PKG.PROB_BOX_FG = 0.0005;
PKG.BOXTVWEIGHT = 0.5;
PKG.DO_BOX_WARP_WEIGHT = true;
PKG.UNARY_LAMBDA = 0.00;
PKG.DO_CONS_PERTURB = true;
PKG.DO_CONS_NOW_PERTURB = true;
PKG.DO_CONSTRAINT_WARP_WEIGHT = true;
PKG.WARP_SAFESPEEDSQUARED = 50.0;
PKG.CONSTRAINT_MEMORY_FACTOR = 0.20;
PKG.OCC_INSIDE_UNITY_MEMORY_FACTOR = 0.50;
PKG.UNITY_MEMORY_FACTOR = 0.20;
PKG.TV_MEMORY_FACTOR = 0.20;
PKG.FG_MEMORY_FACTOR = 0.20;
PKG.DO_UNITY_WARP_WEIGHT = true;
PKG.UNITYHELP = true;
PKG.PROB_UNITY = 0.5;

PKG.CAUSAL = true;
PKG.DO_FORBACKCAUSAL = false; % if true, for noncausal version
PKG.TEST = false; % if true, doesn't process first and last frame
PKG.VIS = 499; % constrols amount of visualization is displayed

if strcmp(model, 'fullboxfig');
  PKG.DO_FORBACKCAUSAL = false;
  outpath = fullfile(plotbase, 'cvos-r1.0-fullboxfig');
elseif strcmp(model, 'fullnopcbf'); % full model but no flow extrapolation
  PKG.DO_CROSSBILATERALFILTERFLOW = false;
  PKG.DO_CONS_PERTURB = false;
  PKG.DO_CONS_NOW_PERTURB = false;
  outpath = fullfile(plotbase, 'wcr-cvos-r1.0-fullnopcbf');
elseif strcmp(model, 'ayv');
  PKG.CAUSAL = true;
  PKG.WEIGHTSHELP = true;
  PKG.UNITYHELP = false;
  PKG.BOXHELP = false;
  PKG.PROB_FG = 0;
  PKG.DO_FORBACKCAUSAL = false;
  PKG.DO_CONS_PERTURB = false;
  PKG.DO_CONS_NOW_PERTURB = false;
  outpath = fullfile(plotbase, 'wcr-cvos-r1.0-ayv');
elseif strcmp(model, 'wfg');
  PKG.CAUSAL = true;
  PKG.WEIGHTSHELP = false;
  PKG.UNITYHELP = false;
  PKG.BOXHELP = false;
  PKG.PROB_FG = 0.0005;
  PKG.DO_FORBACKCAUSAL = false;
  PKG.DO_CONS_PERTURB = false;
  PKG.DO_CONS_NOW_PERTURB = false;
  outpath = fullfile(plotbase, 'wcr-cvos-r1.0-wfg');
elseif strcmp(model, 'fxfperturb');
  PKG.CAUSAL = false;
  PKG.WEIGHTSHELP = false;
  PKG.UNITYHELP = false;
  PKG.BOXHELP = false;
  PKG.PROB_FG = 0;
  PKG.DO_CONS_PERTURB = false;
  PKG.DO_CONS_NOW_PERTURB = true;
  PKG.DO_FORBACKCAUSAL = false;
  outpath = fullfile(plotbase, 'cvos-r1.0-fxfperturb');
end
PKG.outpath = outpath;
try
  if ~exist(outpath, 'dir'); createRequiredFolders(outpath); end;
catch e
  display(e); keyboard; exit;
end

%-----------------------------------------------------------------------------
% run the framework
%-----------------------------------------------------------------------------
for sid = seqs;
  s = num2str(sid{1});
  fprintf('========= running seq: %s ==============\n', s);
  cvos(sid{1}, PKG);
  fprintf('========= finished seq: %s =============\n', s);
end
end
