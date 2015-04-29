function [problem, params] = ttaobdDataMoseg(DATA, paramsIn)

outpath = '/pad_local/btay/projects/cdo/';

%----------------------------------------------------------------------------%
% Sequence Specific 
%----------------------------------------------------------------------------%
switch DATA;
case 1001; id = 10; WB = 07; WF = 07; seq = 'cars1';
case 1002; id = 20; WB = 07; WF = 07; seq = 'cars2';
case 1003; id = 10; WB = 07; WF = 07; seq = 'cars3';
case 1004; id = 20; WB = 07; WF = 07; seq = 'cars4';
case 1005; id = 20; WB = 07; WF = 07; seq = 'cars5';
case 1006; id = 20; WB = 07; WF = 07; seq = 'cars6';
case 1007; id = 20; WB = 07; WF = 07; seq = 'cars7';
case 1008; id = 20; WB = 04; WF = 04; seq = 'cars8';
case 1009; id = 20; WB = 07; WF = 07; seq = 'cars9';
case 1010; id = 20; WB = 07; WF = 07; seq = 'cars10';
case 1011; id = 10; WB = 07; WF = 07; seq = 'people1';
case 1012; id = 10; WB = 07; WF = 07; seq = 'people2';
case 1013; id = 10; WB = 04; WF = 04; seq = 'cars8';
% TODO: tennis, marple
end

switch DATA;
case 1001; iid = 10; WWB = 1; WWF = 1; seq = 'cars1';
case 1002; iid = 20; WWB = 1; WWF = 1; seq = 'cars2';
case 1003; iid = 10; WWB = 1; WWF = 1; seq = 'cars3';
case 1004; iid = 20; WWB = 1; WWF = 1; seq = 'cars4';
case 1005; iid = 20; WWB = 1; WWF = 1; seq = 'cars5';
case 1006; iid = 20; WWB = 1; WWF = 1; seq = 'cars6';
case 1007; iid = 20; WWB = 1; WWF = 1; seq = 'cars7';
case 1008; iid = 20; WWB = 1; WWF = 1; seq = 'cars8';
case 1009; iid = 20; WWB = 1; WWF = 1; seq = 'cars9';
case 1010; iid = 20; WWB = 1; WWF = 1; seq = 'cars10';
case 1011; iid = 10; WWB = 1; WWF = 1; seq = 'people1';
case 1012; iid = 10; WWB = 1; WWF = 1; seq = 'people2';
case 1013; iid = 10; WWB = 1; WWF = 1; seq = 'cars8';
% TODO: tennis, marple

% case 1020 to 1078: total of 59 sequences
case 1020; iid = 20; WWB = 1; WWF = 1; seq = 'bear1';
case 1021; iid = 10; WWB = 1; WWF = 1; seq = 'bear2';
case 1022; iid = 20; WWB = 1; WWF = 1; seq = 'camel1';
case 1023; iid = 10; WWB = 1; WWF = 1; seq = 'cars1';
case 1024; iid = 20; WWB = 1; WWF = 1; seq = 'cars2';
case 1025; iid = 20; WWB = 1; WWF = 1; seq = 'cars3';
case 1026; iid = 20; WWB = 1; WWF = 1; seq = 'cars4';
case 1027; iid = 20; WWB = 1; WWF = 1; seq = 'cars5';
case 1028; iid = 20; WWB = 1; WWF = 1; seq = 'cars6';
case 1029; iid = 20; WWB = 1; WWF = 1; seq = 'cars7';
case 1030; iid = 20; WWB = 1; WWF = 1; seq = 'cars8';
case 1031; iid = 10; WWB = 1; WWF = 1; seq = 'cars9';
case 1032; iid = 20; WWB = 1; WWF = 1; seq = 'cars10';
case 1033; iid = 10; WWB = 1; WWF = 1; seq = 'cats1';
case 1034; iid = 20; WWB = 1; WWF = 1; seq = 'cats2';
case 1035; iid = 20; WWB = 1; WWF = 1; seq = 'cats3';
case 1036; iid = 20; WWB = 1; WWF = 1; seq = 'cats4';
case 1037; iid = 20; WWB = 1; WWF = 1; seq = 'cats5';
case 1038; iid = 20; WWB = 1; WWF = 1; seq = 'cats6';
case 1039; iid = 20; WWB = 1; WWF = 1; seq = 'cats7';
case 1040; iid = 20; WWB = 1; WWF = 1; seq = 'dogs1';
case 1041; iid = 10; WWB = 1; WWF = 1; seq = 'dogs2';
case 1042; iid = 20; WWB = 1; WWF = 1; seq = 'ducks1';
case 1043; iid = 10; WWB = 1; WWF = 1; seq = 'farm1';
case 1044; iid = 20; WWB = 1; WWF = 1; seq = 'giraffes1';
case 1045; iid = 20; WWB = 1; WWF = 1; seq = 'goats1';
case 1046; iid = 20; WWB = 1; WWF = 1; seq = 'horses1';
case 1047; iid = 20; WWB = 1; WWF = 1; seq = 'horses2';
case 1048; iid = 20; WWB = 1; WWF = 1; seq = 'horses3';
case 1049; iid = 20; WWB = 1; WWF = 1; seq = 'horses4';
case 1050; iid = 20; WWB = 1; WWF = 1; seq = 'horses5';
case 1051; iid = 10; WWB = 1; WWF = 1; seq = 'horses6';
case 1052; iid = 20; WWB = 1; WWF = 1; seq = 'lion1';
case 1053; iid = 10; WWB = 1; WWF = 1; seq = 'lion2';
case 1054; iid = 20; WWB = 1; WWF = 1; seq = 'marple1';
case 1055; iid = 20; WWB = 1; WWF = 1; seq = 'marple2';
case 1056; iid = 20; WWB = 1; WWF = 1; seq = 'marple3';
case 1057; iid = 20; WWB = 1; WWF = 1; seq = 'marple4';
case 1058; iid = 20; WWB = 1; WWF = 1; seq = 'marple5';
case 1059; iid = 20; WWB = 1; WWF = 1; seq = 'marple6';
case 1060; iid = 20; WWB = 1; WWF = 1; seq = 'marple7';
case 1061; iid = 10; WWB = 1; WWF = 1; seq = 'marple8';
case 1062; iid = 20; WWB = 1; WWF = 1; seq = 'marple9';
case 1063; iid = 10; WWB = 1; WWF = 1; seq = 'marple10';
case 1064; iid = 20; WWB = 1; WWF = 1; seq = 'marple11';
case 1065; iid = 20; WWB = 1; WWF = 1; seq = 'marple12';
case 1066; iid = 20; WWB = 1; WWF = 1; seq = 'marple13';
case 1067; iid = 20; WWB = 1; WWF = 1; seq = 'meerkats1';
case 1068; iid = 20; WWB = 1; WWF = 1; seq = 'people1';
case 1069; iid = 20; WWB = 1; WWF = 1; seq = 'people2';
case 1070; iid = 20; WWB = 1; WWF = 1; seq = 'people3';
case 1071; iid = 10; WWB = 1; WWF = 1; seq = 'people4';
case 1072; iid = 20; WWB = 1; WWF = 1; seq = 'people5';
case 1073; iid = 10; WWB = 1; WWF = 1; seq = 'rabbits1';
case 1074; iid = 20; WWB = 1; WWF = 1; seq = 'rabbits2';
case 1075; iid = 20; WWB = 1; WWF = 1; seq = 'rabbits3';
case 1076; iid = 20; WWB = 1; WWF = 1; seq = 'rabbits4';
case 1077; iid = 20; WWB = 1; WWF = 1; seq = 'rabbits5';
case 1078; iid = 20; WWB = 1; WWF = 1; seq = 'tennis';
end

if ~exist('id', 'var'); id = iid; end;
if ~exist('WB', 'var'); WB = 1; end;
if ~exist('WF', 'var'); WF = 10; end;

%----------------------------------------------------------------------------%
% data path
%----------------------------------------------------------------------------%
dataset = 'moseg'; ext = 'png';
dpath    = sprintf('/plot/vasiliy/CVPR15/data/moseg/%s', seq);
% uvpath   = sprintf('/plot/vasiliy/CVPR15/data/moseg/%s/flow_g4v', seq);
uvpath   = sprintf('/plot/vasiliy/CVPR15/data/moseg/%s/flow', seq);

nameStr = '%s_%03d';
imgNameStr = ['%s_%02d.ppm.' ext];
uvNameStr  = ['%s_%03d.mat'];

expName = 'moseg';
if exist('paramsIn', 'var'); expName = paramsIn.expName; end;

%----------------------------------------------------------------------------%
% dataset parameters
%----------------------------------------------------------------------------%
params.LAMBDA = 2.0; % 0.5; % 0.8, 0.1; - %ayvaci
params.TAU = 2;
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
