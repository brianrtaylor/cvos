% wrapper that loads current image, flows and occlusions (and uses
% neighboring frames to do this)
% USAGE: [I0, uvf, uvb, uvf_rev, uvb_rev, eb, ef] = load_data(sequence, index)
%
% @param: LOADFLOW : if set to 1, load flow, else, compute it on the spot with TUgraz flow
function [I0, uvf, uvb, uvf_rev, uvb_rev, eb, ef, Ib, If] = load_data( ...
  sequence, index, which, LOADFLOW)

if nargin==2; which = 'all'; end
[~, flow_path, img_path, ~, extension] = dataPaths(sequence);

%-------------------------------------------------------------------------
% load image
%-------------------------------------------------------------------------
img_files = dir([img_path '*.' extension]);
if strcmpi(which,'all')
  I0 = imread( [img_path, '/', img_files(index).name] );
  if size(I0,3)==1, I0 = repmat(I0,[1 1 3]); end
elseif strcmpi(which,'noprev')
  I0 = []; % can be obtained with 'I0 = If'
end

if strcmpi(which,'all')
  Ib = imread( [img_path, '/', img_files(index-1).name] );
  if size(Ib,3)==1, Ib = repmat(Ib,[1 1 3]); end
elseif strcmpi(which,'noprev')
  Ib = []; % can be obtained with 'Ib = I0'
end

If = imread( [img_path, '/', img_files(index+1).name] );

% turn grayscale to RGB if need be
if size(If,3)==1; If = repmat(If,[1 1 3]); end

%-------------------------------------------------------------------------
% load current flow
%-------------------------------------------------------------------------
if LOADFLOW;
  flow_files = dir([flow_path '*.mat']);
  load([flow_path, '/', flow_files(index).name]);
  % load previous and next flows:
  env = load([flow_path, '/', flow_files(index+1).name]);
  uvf_rev = env.uvb;
  
  if strcmpi(which,'all')
    env = load([flow_path, '/', flow_files(index-1).name]);
    uvb_rev = env.uvf;
  elseif strcmpi(which,'noprev')
    uvb_rev = []; % can be obtained with: 'uvb_rev = uvf'
  end
else % compute flow

  % TODO: figure out later, since we have sun or ayvaci flow currently

end

%-------------------------------------------------------------------------
% occlusions
%-------------------------------------------------------------------------
if ~exist('eb', 'var');
  eb = zeros( [size(uvb,1), size(uvb,2)] );
end
if ~exist('ef', 'var')
  ef = zeros( [size(uvf,1), size(uvf,2)] );
end
end
