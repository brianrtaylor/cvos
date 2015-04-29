function [problem, params] = ttaobdDataVasiliy(DATA, paramsIn)

outpath = '/pad_local/btay/projects/cdo/';

%----------------------------------------------------------------------------%
% Sequence Specific 
%----------------------------------------------------------------------------%
switch DATA; % first part kinda useless
case 7001; id = 10; WB = 09; WF = 33; seq = 'frey';
case 7002; id = 10; WB = 09; WF = 08; seq = 'car1_stabilized';
case 7003; id = 10; WB = 09; WF = 30; seq = 'car2_stabilized';
case 7004; id = 10; WB = 09; WF = 1486; seq = 'stefano_ski';
end

switch DATA; % first part kinda useless
case 7001; iid = 10; WWB = 1; WWF = 1; seq = 'frey';
case 7002; iid = 10; WWB = 1; WWF = 1; seq = 'car1_stabilized';
case 7003; iid = 10; WWB = 1; WWF = 1; seq = 'car2_stabilized';
case 7004; iid = 10; WWB = 1; WWF = 1; seq = 'stefano_ski';
end

%----------------------------------------------------------------------------%
% data path
%----------------------------------------------------------------------------%
dataset = 'vasiliyetc'; ext = 'png';
dpath    = sprintf('/plot/vasiliy/CVPR15/data/%s', seq);
uvpath   = sprintf('/plot/vasiliy/CVPR15/data/%s/flow', seq);

nameStr = '%s_%03d';
imgNameStr = ['%s_%03d.' ext];
uvNameStr  = ['%s_%03d.mat'];

expName = 'vasiliyetc';
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
