function [problem, params] = ttaobdDataEtc(DATA, paramsIn)

% outpath = '/pad_local/btay/projects/cdo/';
outpath = '/plot/btay/data/yoseg/';

%----------------------------------------------------------------------------%
% Sequence Specific 
%----------------------------------------------------------------------------%
switch DATA; % first part kinda useless
case 2001; id = 10; WB = 09; WF = 33; seq = 'frey';
case 2002; id = 10; WB = 09; WF = 08; seq = 'car1_stabilized';
case 2003; id = 10; WB = 09; WF = 30; seq = 'car2_stabilized';
case 2004; id = 10; WB = 09; WF = 1486; seq = 'stefano_ski';
case 2005; id = 3; WB = 09; WF = 1486; seq = 'cerealbox';

case 2010; id = 2; WB = 01; WF = 01; seq = 'marple12_10_5';
end

switch DATA; % first part kinda useless
case 2001; iid = 10; WWB = 1; WWF = 1; seq = 'frey';
case 2002; iid = 10; WWB = 1; WWF = 1; seq = 'car1_stabilized';
case 2003; iid = 10; WWB = 1; WWF = 1; seq = 'car2_stabilized';
case 2004; iid = 10; WWB = 1; WWF = 1; seq = 'stefano_ski';
case 2005; iid = 3; WWB = 1; WWF = 1; seq = 'cerealbox';

case 2010; iid = 2; WWB = 01; WWF = 01; seq = 'marple12_10_5';
end

%----------------------------------------------------------------------------%
% data path
%----------------------------------------------------------------------------%
dataset = 'etc'; ext = 'png';
% dpath    = sprintf('/plot/vasiliy/CVPR15/data/%s', seq);
% uvpath   = sprintf('/plot/vasiliy/CVPR15/data/%s/flow', seq);
dpath    = sprintf('/plot/btay/data/yoseg/%s', seq);
uvpath   = sprintf('/plot/btay/data/yoseg/%s/flow', seq);

nameStr = '%s_%03d';
imgNameStr = ['%s_%03d.' ext];
uvNameStr  = ['%s_%03d.mat'];

expName = 'etc';
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
