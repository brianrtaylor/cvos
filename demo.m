% TODO: use better example than cars5 because we lose the right car
% TODO: replace cvos_lite_* with some edits and calls to cvos straight out
% TODO: change all names from tao to cvos
% in frame 1

% 1. download data somewhere local
pth = 'example_demo/';
if ~exist(pth, 'dir'); mkdir(pth); end;
cd(pth);
if ~exist('moseg', 'dir') && ~exist('cars5_demo.zip', 'file');
  unix('wget vision.ucla.edu/~btay/proj/cvos/demo/cars5_demo.zip'); % 204 MB
end
if ~exist('moseg', 'dir');
  unix('unzip cars5_demo.zip');
  % unix('rm cars5_demo.zip');
end
cd('..');

% 2. some setup
% mex_clean
mex_setup

% 3. run cvos_run
% paths to results are editable in cdo_aid/dataPaths.m
cvos_run({'cars5'}, 'full', 'moseg');

% 4. see results in "./results" folder
% TODO: visualize results here
