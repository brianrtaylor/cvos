function [problem, params] = ttaobdDataCube(DATA, paramsIn)
%----------------------------------------------------------------------------%
% Sequence Specific 
%----------------------------------------------------------------------------%

% whole sequence
switch DATA;
case 9001; id = 10; WB = 8; WF = 9; seq = 'cube1';
case 9002; id = 10; WB = 8; WF = 9; seq = 'cube2';
case 9003; id = 10; WB = 8; WF = 9; seq = 'cube3';
case 9004; id = 10; WB = 8; WF = 9; seq = 'cube4';
case 9005; id = 20; WB = 18; WF = 139; seq = 'cube5';
case 9006; id = 20; WB = 18; WF = 129; seq = 'cube6';
case 9007; id = 20; WB = 19; WF = 20; seq = 'cube1a'; % like cube1, moves by 1 pixel, +noise
case 9008; id = 20; WB = 19; WF = 20; seq = 'cube2a'; % like cube2, disparate motion, smaller o--> <-O
case 9009; id = 20; WB = 19; WF = 20; seq = 'cube2b'; % like cube2, disparate motion, smaller o-> <--O
case 9010; id = 20; WB = 19; WF = 20; seq = 'cube2c'; % like cube2, disparate motion, smaller o--> O->
case 9011; id = 20; WB = 19; WF = 20; seq = 'cube2d'; % like cube2, disparate motion, smaller <-o <--O
case 9012; id = 20; WB = 19; WF = 20; seq = 'cube2aa'; % like cube2, disparate motion, smaller o--> <-O (full pass)
case 9013; id = 20; WB = 19; WF = 20; seq = 'cube1b'; % like cube1, smaller
case 9014; id = 20; WB = 19; WF = 20; seq = 'cube1c'; % like cube1, but goes right, stops, & then goes right slower
case 9015; id = 10; WB = 09; WF = 10; seq = 'cube2ab'; % like cube2, but faster
case 9016; id = 10; WB = 09; WF = 10; seq = 'cube7'; % like cube 2, o-> OV o.(down, stationary)
end

% single frame
switch DATA;
case 9001; iid = 10; WWB =  1; WWF = 1; seq = 'cube1';
case 9002; iid = 10; WWB =  8; WWF = 9; seq = 'cube2';
case 9003; iid = 10; WWB =  8; WWF = 9; seq = 'cube3';
case 9004; iid = 10; WWB =  8; WWF = 9; seq = 'cube4';
case 9005; iid = 20; WWB = 18; WWF = 139; seq = 'cube5';
case 9006; iid = 20; WWB = 18; WWF = 129; seq = 'cube6';
case 9007; iid = 20; WWB = 3; WWF = 3; seq = 'cube1a'; % like cube1, moves by 1 pixel, +noise
case 9008; iid = 35; WWB = 1; WWF = 1; seq = 'cube2a'; % like cube2, disparate motion, smaller o--> <-O
case 9009; iid = 20; WWB = 19; WWF = 20; seq = 'cube2b'; % like cube2, disparate motion, smaller o-> <--O
case 9010; iid = 20; WWB = 19; WWF = 20; seq = 'cube2c'; % like cube2, disparate motion, smaller o--> O->
case 9011; iid = 20; WWB = 19; WWF = 20; seq = 'cube2d'; % like cube2, disparate motion, smaller <-o <--O
case 9012; iid = 20; WWB = 19; WWF = 20; seq = 'cube2aa'; % like cube2, disparate motion, smaller o--> <-O (full pass)
case 9013; iid = 20; WWB = 19; WWF = 20; seq = 'cube1b'; % like cube1, smaller
case 9014; iid = 5; WWB = 1; WWF = 1; seq = 'cube1c'; % portion of cube1c moving right at speed 3 pxls
case 9015; iid = 12; WWB = 1; WWF = 1; seq = 'cube2ab'; % like cube2, but faster
case 9016; iid = 4; WWB = 1; WWF = 1; seq = 'cube7'; % 
end
%----------------------------------------------------------------------------%
% Dataset wide
%----------------------------------------------------------------------------%
dataset = 'cube'; ext = 'png';
dpath    = '/pad_local/btay/data/cube/images/';
uvpath   = '/pad_local/btay/data/cube/flow/';
ucm2path = '/pad_local/btay/data/cube/gpb/';

expName = 'cube_exp1';
if exist('paramsIn', 'var'); expName = paramsIn.expName; end;

nameStr = '%s_%03d';
imgNameStr = [nameStr '.' ext];
% uvNameStr  = [nameStr '-uv-%s.mat'];
uvNameStr  = [nameStr '-uv-gt.mat']; % ground truth flow here (for loading cdo)
% intuvNameStr  = [nameStr '-uv-%s-b%02d-f%02d.mat'];
% intuvNameStr  = '%s_int_%s_to_%03d_from_%%03d-uv-%s.mat';
% gpbNameStr = nameStr;
% resNameStr = [nameStr '-b%02d-f%02d_' expName];

gtName = '/pad_local/btay/data/cube/GT/GT.mat';
params.outpath = '/pad_local/btay/projects/taobd/';

%----------------------------------------------------------------------------%
% PARAMETERS: TODO: check these paramteres, shrink them a lot! 
%----------------------------------------------------------------------------%
% recomputing
params.RECOMPUTE = false;
params.RECOMPUTE_INTEGRATED_FLOW = false;
params.RECOMPUTE_INTERFRAME_FLOW = false;
params.RECOMPUTE_GPB = false;
params.RECOMPUTE_UCM2 = false;

params.LAMBDA = 0.2;
params.TAU = 1;
params.OCCPROB = 0.90;
params.MINCONSTRAINTDIST = 2.0;

%----------------------------------------------------------------------------%
% Outputs
%----------------------------------------------------------------------------%
if exist('paramsIn', 'var');
  params = structmerge(params, paramsIn);
end

% build problem for dataset
% problem = v2struct(seq, ext, id, WB, WF, nameStr, imgNameStr, uvNameStr, ...
%   gpbNameStr, resNameStr, gtName, intuvNameStr, dpath, uvpath, ucm2path, ...
%   dataset, expName);
problem = v2struct(seq, ext, id, WB, WF, iid, WWB, WWF, ...
  nameStr, imgNameStr, uvNameStr, gtName, dpath, uvpath, ucm2path, ...
  dataset, expName);
end
