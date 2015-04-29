function [problem, params] = ttaobdDataSegtrack(DATA, paramsIn)

outpath = '/pad_local/btay/projects/cdo/';

%----------------------------------------------------------------------------%
% Sequence Specific 
%----------------------------------------------------------------------------%
switch DATA;
case 4001; id = 10; WB = 07; WF = 07; seq = 'birdfall2';
case 4002; id = 20; WB = 07; WF = 07; seq = 'cheetah';
case 4003; id = 10; WB = 07; WF = 07; seq = 'girl';
case 4004; id = 20; WB = 07; WF = 07; seq = 'monkeydog';
case 4005; id = 20; WB = 07; WF = 07; seq = 'parachute';
case 4006; id = 20; WB = 07; WF = 07; seq = 'penguin';
end

switch DATA;
case 4001; iid = 10; WWB = 1; WWF = 1; seq = 'birdfall2';
case 4002; iid = 20; WWB = 1; WWF = 1; seq = 'cheetah';
case 4003; iid = 10; WWB = 1; WWF = 1; seq = 'girl';
case 4004; iid = 20; WWB = 1; WWF = 1; seq = 'monkeydog';
case 4005; iid = 20; WWB = 1; WWF = 1; seq = 'parachute';
case 4006; iid = 20; WWB = 1; WWF = 1; seq = 'penguin';
end

%----------------------------------------------------------------------------%
% data path
%----------------------------------------------------------------------------%
dataset = 'segtrack'; ext = 'png';
dpath    = sprintf('/plot/vasiliy/CVPR15/data/segtrack/%s', seq);
uvpath   = sprintf('/plot/vasiliy/CVPR15/data/segtrack/%s/flow', seq);

nameStr = '%s_%03d';
imgNameStr = ['%s_%03d.' ext];
uvNameStr  = ['%s_%03d.mat'];

expName = 'segtrack';
if exist('paramsIn', 'var'); expName = paramsIn.expName; end;

%----------------------------------------------------------------------------%
% dataset parameters
%----------------------------------------------------------------------------%
params.LAMBDA = 1.0; % 0.5; % 0.8, 0.1; - %ayvaci
params.TAU = 1e-4;
params.OCCPROB = 0.50;

%----------------------------------------------------------------------------%
% Outputs
%----------------------------------------------------------------------------%
if exist('paramsIn', 'var');
  params = structmerge(params, paramsIn);
end

problem = v2struct(seq, ext, id, WB, WF, iid, WWB, WWF, ...
  nameStr, imgNameStr, uvNameStr, dpath, uvpath, dataset, expName);
end
