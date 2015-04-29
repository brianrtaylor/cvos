function [problem, params] = ttaobdDataCmuOccFull(DATA, paramsIn)

outpath = '/pad_local/btay/projects/taobd/batchrun2/';

%----------------------------------------------------------------------------%
% Sequence Specific 
%----------------------------------------------------------------------------%
switch DATA;
% testing this bench one
case 6001; id = 10; WB = 10; WF = 09; cr = [5 7 311 238]; seq = 'bench'; % goes from image 0 to 19
case 6002; id = 10; WB = 10; WF = 08; cr = [25 6 312 236]; seq = 'car'; % img19 has annoying specularity (ignore)
case 6003; id = 10; WB = 10; WF = 09; cr = [6 3 317 238]; seq = 'chair1'; % img0 to img1 has large motion (keep!)
case 6004; id = 10; WB = 10; WF = 09; cr = [4 7 317 234]; seq = 'cmu_sign';
case 6005; id = 10; WB = 10; WF = 09; cr = [8 7 310 237]; seq = 'coffee_stuff';
case 6006; id = 10; WB = 10; WF = 09; cr = [9 3 317 238]; seq = 'couch_color';
case 6007; id = 05; WB = 05; WF = 04; cr = [4 3 318 238]; seq = 'couch_corner';
case 6008; id = 10; WB = 10; WF = 09; cr = [4 7 317 238]; seq = 'fencepost'; % img0 goes different direction 
case 6009; id = 15; WB = 15; WF = 14; cr = [4 4 318 237]; seq = 'hand2'; % img28 has green line (change WF to 12?)
case 6010; id = 05; WB = 05; WF = 04; cr = [4 3 316 238]; seq = 'hand3';
% case 6011; id = 10; WB = 09; WF = 09; cr = [18 5 313 235]; seq = 'intrepid';
% case 6012; id = 10; WB = 09; WF = 09; cr = [3 13 306 230]; seq = 'intrepid_corner';
% case 6013; id = 10; WB = 09; WF = 09; cr = [4 3 311 238]; seq = 'intrepid_corner2';
case 6011; id = 10; WB = 10; WF = 09; cr = [3 13 306 230]; seq = 'intrepid_corner'; % img0 has red line 
case 6012; id = 10; WB = 10; WF = 09; cr = [4 3 311 238]; seq = 'intrepid_corner2'; % added img0
case 6013; id = 10; WB = 10; WF = 09; cr = [18 5 313 235]; seq = 'intrepid'; % added img0
case 6014; id = 10; WB = 10; WF = 09; cr = [3 6 314 238]; seq = 'linus1'; % added img0 (~no motion)
case 6015; id = 15; WB = 15; WF = 14; cr = [3 3 317 236]; seq = 'mugs'; % img3 has blue line
case 6016; id = 10; WB = 10; WF = 09; cr = [8 3 312 233]; seq = 'mugs2'; % added img0 (good motion)
case 6017; id = 10; WB = 10; WF = 09; cr = [8 4 312 237]; seq = 'post';
case 6018; id = 10; WB = 10; WF = 09; cr = [3 5 316 229]; seq = 'rocking_horse'; % added img0 (~no motion)
case 6019; id = 04; WB = 04; WF = 03; cr = [3 3 316 208]; seq = 'squirrel2'; % added img0 (good motion)
case 6020; id = 05; WB = 05; WF = 04; cr = [3 4 312 165]; seq = 'squirrel3'; % img0 good motion, kinda blurry
case 6021; id = 04; WB = 04; WF = 03; cr = [6 4 317 159]; seq = 'squirrel4'; % added img0 (good motion)
case 6022; id = 10; WB = 10; WF = 09; cr = [8 5 317 238]; seq = 'staplers'; % added img0
case 6023; id = 10; WB = 10; WF = 09; cr = [20 5 303 220]; seq = 'trash'; % added img0 (good motion)
case 6024; id = 10; WB = 10; WF = 09; cr = [8 3 308 238]; seq = 'trash_can'; % added img0 (good motion)
case 6025; id = 10; WB = 10; WF = 09; cr = [7 6 314 238]; seq = 'tree'; % added img0
case 6026; id = 05; WB = 05; WF = 04; cr = [3 3 317 238]; seq = 'walking_legs'; % added img0 (huge motion)
case 6027; id = 10; WB = 10; WF = 09; cr = [5 3 315 233]; seq = 'wooden_man'; % added img0 (good motion)
% case 6028; id = 05; WB = 05; WF = 04; cr = [6 3 315 234]; seq = 'zoe1'; img9 is same as img0, img1 (ignore it)
case 6028; id = 05; WB = 05; WF = 03; cr = [6 3 315 234]; seq = 'zoe1'; % added img0, img9 == img0, img1 (ignore)
case 6029; id = 10; WB = 10; WF = 09; cr = [14 19 301 237]; seq = 'zoe2'; % added img0 (good motion)
case 6030; id = 06; WB = 06; WF = 05; cr = [4 4 317 237]; seq = 'zoe3'; % added img0 (huge motion)
end

switch DATA;
% testing this bench one
case 6001; iid = 10; WWB = 10; WWF = 09; seq = 'bench'; % goes from image 0 to 19
case 6002; iid = 10; WWB = 10; WWF = 08; seq = 'car'; % img19 has annoying specularity (ignore)
case 6003; iid = 10; WWB = 10; WWF = 09; seq = 'chair1'; % img0 to img1 has large motion (keep!)
case 6004; iid = 10; WWB = 10; WWF = 09; seq = 'cmu_sign';
case 6005; iid = 10; WWB = 10; WWF = 09; seq = 'coffee_stuff';
case 6006; iid = 10; WWB = 10; WWF = 09; seq = 'couch_color';
case 6007; iid = 05; WWB = 05; WWF = 04; seq = 'couch_corner';
case 6008; iid = 10; WWB = 10; WWF = 09; seq = 'fencepost'; % img0 goes different direction 
case 6009; iid = 15; WWB = 15; WWF = 14; seq = 'hand2'; % img28 has green line (change WF to 12?)
case 6010; iid = 05; WWB = 05; WWF = 04; seq = 'hand3';
% case 6011; iid = 10; WWB = 09; WWF = 09; seq = 'intrepid';
% case 6012; iid = 10; WWB = 09; WWF = 09; seq = 'intrepid_corner';
% case 6013; iid = 10; WWB = 09; WWF = 09; seq = 'intrepid_corner2';
case 6011; iid = 10; WWB = 10; WWF = 09; seq = 'intrepid_corner'; % img0 has red line 
case 6012; iid = 10; WWB = 10; WWF = 09; seq = 'intrepid_corner2'; % added img0
case 6013; iid = 10; WWB = 10; WWF = 09; seq = 'intrepid'; % added img0
case 6014; iid = 10; WWB = 10; WWF = 09; seq = 'linus1'; % added img0 (~no motion)
case 6015; iid = 15; WWB = 15; WWF = 14; seq = 'mugs'; % img3 has blue line
case 6016; iid = 10; WWB = 10; WWF = 09; seq = 'mugs2'; % added img0 (good motion)
case 6017; iid = 10; WWB = 10; WWF = 09; seq = 'post';
case 6018; iid = 10; WWB = 10; WWF = 09; seq = 'rocking_horse'; % added img0 (~no motion)
case 6019; iid = 04; WWB = 04; WWF = 03; seq = 'squirrel2'; % added img0 (good motion)
case 6020; iid = 05; WWB = 05; WWF = 04; seq = 'squirrel3'; % img0 good motion, kinda blurry
case 6021; iid = 04; WWB = 04; WWF = 03; seq = 'squirrel4'; % added img0 (good motion)
case 6022; iid = 10; WWB = 10; WWF = 09; seq = 'staplers'; % added img0
case 6023; iid = 10; WWB = 10; WWF = 09; seq = 'trash'; % added img0 (good motion)
case 6024; iid = 10; WWB = 10; WWF = 09; seq = 'trash_can'; % added img0 (good motion)
case 6025; iid = 10; WWB = 10; WWF = 09; seq = 'tree'; % added img0
case 6026; iid = 05; WWB = 05; WWF = 04; seq = 'walking_legs'; % added img0 (huge motion)
case 6027; iid = 10; WWB = 10; WWF = 09; seq = 'wooden_man'; % added img0 (good motion)
% case 6028; iid = 05; WWB = 05; WWF = 04; seq = 'zoe1'; img9 is same as img0, img1 (ignore it)
case 6028; iid = 05; WWB = 05; WWF = 03; seq = 'zoe1'; % added img0, img9 == img0, img1 (ignore)
case 6029; iid = 10; WWB = 10; WWF = 09; seq = 'zoe2'; % added img0 (good motion)
case 6030; iid = 06; WWB = 06; WWF = 05; seq = 'zoe3'; % added img0 (huge motion)
end

% alper say tests sequences are as follows:
% bench, car2. chair1. coffee stuff, couch color, couch corner, fencepost, hand3, intrepid, post, rocking horse, squirrel4, trash can, tree, walking legs, zoe1

%----------------------------------------------------------------------------%
% data path
%----------------------------------------------------------------------------%
dataset = 'cmufull'; ext = 'tif';
dpath = sprintf('/pad_local/btay/data/cmu-occlusion/images/raw/%s/', seq);
uvpath  = sprintf('/pad_local/btay/data/cmu-occlusion/flow/raw/%s/', seq);

nameStr = '%s_%03d';
imgNameStr = [nameStr '.' ext];
uvNameStr  = [nameStr '.mat'];

expName = 'play';
if exist('paramsIn', 'var'); expName = paramsIn.expName; end;

%----------------------------------------------------------------------------%
% dataset parameters
%----------------------------------------------------------------------------%
params.LAMBDA = 1.0; % 20131010
params.TAU = 1.0; % 20131010
params.OCCPROB = 0.1;

%----------------------------------------------------------------------------%
% Outputs
%----------------------------------------------------------------------------%
if exist('paramsIn', 'var');
  params = structmerge(params, paramsIn);
end

problem = v2struct(seq, ext, id, WB, WF, iid, WWB, WWF, ...
  nameStr, imgNameStr, uvNameStr, dpath, uvpath, dataset, expName);
end
