%-----------------------------------------------------------------------------
% dataPaths
%
% contains default paths to the datasets we worked with. Please change 
% variable BASE_PATH to the appropriate one in your directory structure
%
% If you want to use a new sequence, the easiest way is to format the
% data as we have here by placing the image and flow files into two 
% separate directories and providings their paths and extensions to this file
%
% @return: sequence_name: sequence prefix or name
% @return: flow_path: path to the flow files
% @return: img_path: path to the image files
% @return: truth_path: path to the groundtruth files
% @return: extension: image file extension
% @return: flowtype: if you have different methods for computing optical flow 
%   and wish to compare their effects on the performance of this algorithm
% @param: seq: sequence name (needs to be unique across all datasets you use)
%-----------------------------------------------------------------------------
function [sequence_name, flow_path, img_path, truth_path, extension, ...
  flowtype] = dataPaths(seq)

BASE_PATH = 'example_demo/';

sequence_name = seq;
flow_path = '';
img_path = '';
truth_path = '';
extension = '';
flowtype = '';

data = set_data_struct(BASE_PATH);
h = cat(1, data(:).hash );
idx = find( string2hash(seq) == h);
if isempty(idx)
  fprintf('Did not find "%s" in the dataset. Is the name mispelled?\n', seq);
  return;
end

sequence_name = data(idx).sequence_name;
flow_path = data(idx).flow_path;
img_path = data(idx).img_path;
truth_path = data(idx).truth_path;
extension = data(idx).extension;
flowtype = data(idx).flowtype;
end

%-----------------------------------------------------------------------------
% helper function with all the info
%-----------------------------------------------------------------------------
function data = set_data_struct(BASE_PATH)
k=1;
data(k).truth_path = '';
% -------------------------------------------------------------------------
%                           MOSEG dataset
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars3';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars4';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars5';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars6';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars7';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars8';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars9';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cars10';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'people1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'people2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple3';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple4';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple5';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple6';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple6';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple7';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple8';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple9';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple10';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple11';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple12';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'marple13';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'tennis';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'camel1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'dogs1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'dogs2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'farm1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'goats1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'horses1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'horses2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'horses3';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'horses4';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'horses5';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'horses6';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'lion1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'lion2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'people3';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'people4';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'people5';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'giraffes1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'rabbits1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'rabbits2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'rabbits3';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'rabbits4';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'rabbits5';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'bear1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'bear2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cats1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cats2';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cats3';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cats4';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cats5';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cats6';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'cats7';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'ducks1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
% -------------------------------------------------------------------------
data(k).sequence_name = 'meerkats1';
data(k).flow_path = [BASE_PATH, '/moseg/', data(k).sequence_name, '/flow/'];
data(k).img_path = [BASE_PATH,'/moseg/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                  Berkeley Video Segmentation Dataset
% -------------------------------------------------------------------------
data(k).sequence_name = 'hockey_goals';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'hummingbird';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'belly_dancing';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'rock_climbing';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'kim_yu_na';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'animal_chase';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'rock_climbing2';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'buck';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'humpback';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'yosemite';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'penguins';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'panda';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'snow_leopards';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'buffalos';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'frozen_lake';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'planet_earth_2';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'horse_gate';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'salsa';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'nordic_skiing';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'kangaroo_fighting';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'snowboarding';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'fisheye';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'street_food';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'palm_tree';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'snowboarding_crashes';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'hockey';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'monkeys_behind_fence';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'fish_underwater';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'swimming';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'beach_volleyball';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'arctic_kayak';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'birds_of_paradise';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'shark_attack';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'gokart';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'ballet';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'capoeira';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'nba_commercial';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'samba_kids';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'riverboat';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'bicycle_race';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'new_york';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'white_tiger';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'octopus';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'harley_davidson';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'jungle_cat';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'freight_train';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'baseball';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'koala';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'up_trailer';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'beyonce';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'airplane';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'pepsis_wasps';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'planet_earth_1';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'chameleons';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'slow_polo';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'up_dug';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'panda_cub';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'sled_dog_race';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'angkor_wat';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------
data(k).sequence_name = 'vw_commercial';
data(k).flow_path = [BASE_PATH, '/bvds/test/', data(k).sequence_name, '/flow_small/'];
data(k).img_path = [BASE_PATH,'/bvds/test/', data(k).sequence_name, '/'];
data(k).extension = 'png';
data(k).flowtype = 'sun';
k=k+1;
%--------------------------------------------------------------------

% compute hashes for strings
for ii=1:length(data)
   data(ii).hash = string2hash( data(ii).sequence_name ); 
end
end
