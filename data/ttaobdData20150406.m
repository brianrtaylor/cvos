function [problem, params] = ttaobdData(DATA, paramsIn)

if ~exist('paramsIn', 'var'); paramsIn = struct; paramsIn.expName = 'test'; end;
dataset = 'anon';
expName = 'anon';
params = struct;

%----------------------------------------------------------------------------%
% Some parameters
%----------------------------------------------------------------------------%
params.FILTERDIVUV_STEPUV = false;
params.FILTERDIVUV_INTUV = false; % does nothing for now
params.FILTERDIVUV_OCCCUES = false;
params.REFINED_FLOW = false;
params.AFFINE_FLOW = true;
params.FLOWOPTS = ''; % '-decomp'
params.FLOWMETHOD = ['ayvaci-gen' params.FLOWOPTS]; % 'ayvaci'; % old parameters ayvaci don't forget
% params.FLOWMETHOD = ['deqing' params.FLOWOPTS]; % 'ayvaci';
% params.FLOWMETHOD = ['brox2alp' params.FLOWOPTS]; % 'ayvaci';
% params.FLOWMETHOD = ['brox' params.FLOWOPTS]; % 'ayvaci';
% params.FLOWMETHOD = ['ayvaci-kr' params.FLOWOPTS]; % 'ayvaci';

% occ detection
params.DETECT_OCC_METHOD = 'THRESHOLD_LAB_RES';
% params.DETECT_OCC_METHOD = 'THRESHOLD_RGB_RES';
% params.DETECT_OCC_METHOD = 'THRESHOLD_RGB_RES_ADAPTIVE';
params.OCC_LAB_THRESH = 15;
params.OCC_RGB_THRESH = 15; % stupidly used for LAB threshold
% params.OCC_RGB_THRESH = 50;

% flow confidence
params.SIGMOID_LOW_SLOPE = 10;
params.SIGMOID_HIGH_SLOPE = 10;
params.SIGMOID_LOW_THRESH = 2;
params.SIGMOID_HIGH_THRESH = 15;
% params.SEGTREE_THRESHOLDS = [0.02];
params.SEGTREE_THRESHOLDS = [0.01];
% params.SEGTREE_THRESHOLDS = [0.001];
% params.LAMBDA = 0.5; % 0.8, 0.1; - %ayvaci
% params.LAMBDA = 0.1;
% params.TAU = 3;
% params.AGGMETHOD = 'alper'; %btay
% params.CUEMETHOD = 'alper'; %btay
% params.SCOREMETHOD = 'alper'; %btay
params.AGGMETHOD = 'btay-best';
params.CUEMETHOD = 'btay';
params.SCOREMETHOD = 'btay';

params.NEGDIV_COEFF = 0.0;
params.POSDIV_COEFF = 1.0;

params.RECOMPUTE = false;
params.RECOMPUTE_INTEGRATED_FLOW = true;
params.RECOMPUTE_INTERFRAME_FLOW = false;
params.RECOMPUTE_GPB = false;
params.RECOMPUTE_UCM2 = false;

params.MINCONSTRAINTDIST = 2.0;

% visuals
params.VIS_FLOW = false;
params.VIS_OCC = false;
params.VIS_DIVUV = false;
params.VIS_CUES = false;
params.SELFEVAL = false;

params.NEGDIVCUE = 0;



%===============================================================================
%                                                                      DATA PATH
%===============================================================================
switch DATA;

case 0; % cube sequence
dataset = 'cube';
seq = 'cube1'; ext = 'png'; id = 10; WF = 10; WB = 9;
expName = 'cube_exp1';

if exist('paramsIn', 'var') && isfield(paramsIn, 'expName');
  expName = paramsIn.expName;
end

dpath    = '/pad/btay/data/cube/images/';
uvpath   = '/pad/btay/data/cube/flow/';
ucm2path = '/pad/btay/data/cube/gpb/';

nameStr = '%s_%03d';
imgNameStr = [nameStr '.' ext];
uvNameStr  = [nameStr '-uv-%s.mat'];
intuvNameStr  = [nameStr '-uv-%s-b%02d-f%02d.mat'];
gpbNameStr = nameStr;
resNameStr = [nameStr '-b%02d-f%02d_' expName];

gtName = '/pad/btay/data/cube/GT/GT.mat';

params.outpath = '/pad/btay/projects/taobd/';

params.SIGMOID_LOW_SLOPE = 5;
params.SIGMOID_HIGH_SLOPE = 5;
params.SIGMOID_LOW_THRESH = 2;
params.SIGMOID_HIGH_THRESH = 15;
params.SEGTREE_THRESHOLDS = [0.02];
% params.LAMBDA = 0.2;
% params.TAU = 1;
params.SELFEVAL = false;

case 1; % cube2 sequence
dataset = 'cube';
% seq = 'cube2'; ext = 'png'; id = 10; WF = 10; WB = 9;
seq = 'cube4'; ext = 'png'; id = 10; WF = 9; WB = 9;
expName = 'cube_exp4';

if exist('paramsIn', 'var') && isfield(paramsIn, 'expName');
  expName = paramsIn.expName;
end

dpath    = '/pad/btay/data/cube/images/';
uvpath   = '/pad/btay/data/cube/flow/';
ucm2path = '/pad/btay/data/cube/gpb/';

nameStr = '%s_%03d';
imgNameStr = [nameStr '.' ext];
uvNameStr  = [nameStr '-uv-%s.mat'];
intuvNameStr  = [nameStr '-uv-%s-b%02d-f%02d.mat'];
gpbNameStr = nameStr;
resNameStr = [nameStr '-b%02d-f%02d_' expName];

gtName = '/pad/btay/data/cube/GT/GT.mat';

params.outpath = '/pad/btay/projects/taobd/';

params.FLOWMETHOD = 'gt';

params.SIGMOID_LOW_SLOPE = 5;
params.SIGMOID_HIGH_SLOPE = 10;
params.SIGMOID_LOW_THRESH = 2;
params.SIGMOID_HIGH_THRESH = 15;
params.SEGTREE_THRESHOLDS = [0.02];
% params.LAMBDA = 0.2;
% params.TAU = 1;
params.SELFEVAL = false;

case 9; % Soccer sequnce
dataset = 'soccer';
seq = 'soccer'; ext = 'png'; id = 162; WF = 15; WB = 15;
% imgNameStr = '%s-%05d.png'; uvNameStr = '%s-%05d-uv.mat';
% dpath    = '/scratch/ayvaci/vidsegCVPR11/soccer/frames/';
% uvpath   = '/scratch/ayvaci/vidsegCVPR11/soccer/flow/';
% ucm2path = '/scratch/ayvaci/vidsegCVPR11/soccer/ucm2/';
imgNameStr = '%s_%06d.png'; uvNameStr = '%s_%06d-uv.mat';
gpbNameStr = '%s_%06d.pb'; resNameStr = '%s_%06d.mat';
dpath    = '/pad/btay/data/soccer/images/';
uvpath   = '/pad/btay/data/soccer/flow/pm1/';
ucm2path = '/pad/btay/data/soccer/gpb/';

expName = 'soccerExp';
nameStr = '%s_%06d';
imgNameStr = [nameStr '.' ext];
uvNameStr  = [nameStr '-uv-%s.mat'];
intuvNameStr  = [nameStr '-uv-%s-b%02d-f%02d.mat'];
gpbNameStr = nameStr;
resNameStr = [nameStr '-b%02d-f%02d_' expName];
params.outpath = '/pad/btay/projects/taobd/soccerFig/';
% exp1: no magic filter at all cmu

% groundtruth
gtName = sprintf('%s-gt', sprintf(nameStr, seq, id));


params.SIGMOID_LOW_SLOPE = 5;
params.SIGMOID_HIGH_SLOPE = 10;
params.SIGMOID_LOW_THRESH = 2;
params.SIGMOID_HIGH_THRESH = 15;
params.SEGTREE_THRESHOLDS = [0.02];
% params.LAMBDA = 0.2;
% params.TAU = 1;
params.SELFEVAL = false;

end




%===============================================================================
%                                                               SOME NASTY HACKS 
%===============================================================================

% intuvNameStr defined across whole dataset, there are a few others like this too
% should probably conform to this style perhaps
if strcmp(dataset, 'bvds');
  intuvNameStr = '%s_%03d-uv-%02d.mat';
elseif strcmp(dataset, 'camvid');
  intuvNameStr = '%s_%06d-uv-%02d.mat';
elseif strcmp(dataset, 'moseg');
  intuvNameStr = '%s_%06d-uv-%02d.mat';
else
  intuvNameStr = '%s_%06d-uv-%02d.mat';
end


%===============================================================================
%                                                                   OUTPUTS
%===============================================================================

% do something else just in cases
% if it's cmu-occlusion, get it from that file
if (DATA > 1000) && (DATA < 1079);
  [problem, params] = ttaobdDataMoseg(DATA, paramsIn);
elseif (DATA > 2000) && (DATA < 2010);
  % [problem, params] = ttaobdDataMsee2014(DATA, paramsIn);
  [problem, params] = ttaobdDataEtc(DATA, paramsIn);
elseif (DATA > 3000) && (DATA < 3700);
  [problem, params] = ttaobdDataCamvid(DATA, paramsIn);
elseif (DATA > 4000) && (DATA < 4010);
  [problem, params] = ttaobdDataSegtrack(DATA, paramsIn);
elseif (DATA > 5000) && (DATA < 5016);
  [problem, params] = ttaobdDataVideosegmentation(DATA, paramsIn);
elseif (DATA > 6000) && (DATA < 6031);
  [problem, params] = ttaobdDataCmuOccFull(DATA, paramsIn);
elseif (DATA > 7000) && (DATA < 7010);
  [problem, params] = ttaobdDataVasiliy(DATA, paramsIn);
  % [problem, params] = ttaobdDataCmuOcc(DATA, paramsIn);
elseif (DATA > 8000) && (DATA < 8201);
  [problem, params] = ttaobdDataBvds(DATA, paramsIn);
elseif (DATA > 9000) && (DATA < 9020);
  [problem, params] = ttaobdDataCube(DATA, paramsIn);
else
  if ~isfield(params, 'FLOWMETHOD');
    params.FLOWMETHOD = 'ayvaci-gen';
  end
  
  % build problem for dataset
  problem = v2struct(DATA, seq, ext, id, WF, WB, nameStr, imgNameStr, ...
    uvNameStr, gpbNameStr, resNameStr, gtName, intuvNameStr, ...
    dpath, uvpath, ucm2path, dataset, expName);
end

% build params
if exist('paramsIn', 'var');
  params = structmerge(params, paramsIn);
end

end
