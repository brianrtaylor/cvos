%-----------------------------------------------------------------------------
% demo
%
% this code downloads a sequence from the website and then runs the cvos
% framework on that sequence, storing the results in the results folder
%-----------------------------------------------------------------------------

% 1. download data
pth = 'example_demo/';
if ~exist(pth, 'dir'); mkdir(pth); end;
cd(pth);
if ~exist('moseg', 'dir') && ~exist('cars5_demo.zip', 'file');
  unix('wget vision.ucla.edu/~btay/proj/cvos/demo/cars5_demo.zip'); % ~204 MB
end
if ~exist('moseg', 'dir');
  unix('unzip cars5_demo.zip');
  % unix('rm cars5_demo.zip');
end
cd('..');

% 2. some setup
mex_clean
mex_setup

% 3. run cvos_run
% paths to results are editable in cvos_aid/dataPaths.m
cvos_run({'cars5'}, 'full');

% 4. checkout the results in "./results" folder!
