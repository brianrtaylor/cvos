function [problem, params] = ttaobdDataVideosegmentation(DATA, paramsIn)

outpath = '/pad_local/btay/projects/cdo/';

%----------------------------------------------------------------------------%
% Sequence Specific 
%----------------------------------------------------------------------------%
switch DATA; % first part kinda useless
case 5001; id = 10; WB = 09; WF = 20; seq = 'atonement';
case 5002; id = 10; WB = 09; WF = 20; seq = 'coraline';
case 5003; id = 10; WB = 09; WF = 20; seq = 'diving';
case 5004; id = 10; WB = 09; WF = 20; seq = 'earth';
case 5005; id = 10; WB = 09; WF = 20; seq = 'flower_garden';
case 5006; id = 10; WB = 09; WF = 20; seq = 'football';
case 5007; id = 10; WB = 09; WF = 20; seq = 'goodfellas_kitchen';
case 5008; id = 10; WB = 09; WF = 20; seq = 'kimyuna2';
case 5009; id = 10; WB = 09; WF = 20; seq = 'nocountryforoldmen';
case 5010; id = 10; WB = 09; WF = 20; seq = 'paris_compare';
case 5011; id = 10; WB = 09; WF = 20; seq = 'publicenemies1';
case 5012; id = 10; WB = 09; WF = 20; seq = 'publicenemies2';
case 5013; id = 10; WB = 09; WF = 20; seq = 'slomo_surfer';
case 5014; id = 10; WB = 09; WF = 20; seq = 'slumdog1';
case 5015; id = 10; WB = 09; WF = 20; seq = 'waterski';
end

switch DATA; % first part kinda useless
case 5001; iid = 10; WWB = 1; WWF = 1; seq = 'atonement';
case 5002; iid = 10; WWB = 1; WWF = 1; seq = 'coraline';
case 5003; iid = 10; WWB = 1; WWF = 1; seq = 'diving';
case 5004; iid = 10; WWB = 1; WWF = 1; seq = 'earth';
case 5005; iid = 10; WWB = 1; WWF = 1; seq = 'flower_garden';
case 5006; iid = 10; WWB = 1; WWF = 1; seq = 'football';
case 5007; iid = 10; WWB = 1; WWF = 1; seq = 'goodfellas_kitchen';
case 5008; iid = 10; WWB = 1; WWF = 1; seq = 'kimyuna2';
case 5009; iid = 10; WWB = 1; WWF = 1; seq = 'nocountryforoldmen';
case 5010; iid = 10; WWB = 1; WWF = 1; seq = 'paris_compare';
case 5011; iid = 10; WWB = 1; WWF = 1; seq = 'publicenemies1';
case 5012; iid = 10; WWB = 1; WWF = 1; seq = 'publicenemies2';
case 5013; iid = 10; WWB = 1; WWF = 1; seq = 'slomo_surfer';
case 5014; iid = 10; WWB = 1; WWF = 1; seq = 'slumdog1';
case 5015; iid = 10; WWB = 1; WWF = 1; seq = 'waterski';
end

%----------------------------------------------------------------------------%
% data path
%----------------------------------------------------------------------------%
dataset = 'videosegmentation'; ext = 'png';
dpath    = sprintf('/plot/vasiliy/CVPR15/data/videosegmentation/%s', seq);
uvpath   = sprintf('/plot/vasiliy/CVPR15/data/videosegmentation/%s/flow', seq);

nameStr = '%s_%03d';
imgNameStr = ['%s_%03d.' ext];
uvNameStr  = ['%s_%03d.mat'];

expName = 'videosegmentation';
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
