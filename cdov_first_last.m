%----------------------------------------------------------------------------%
% cdov(DATA, PKG)
%
% causal detachable objects
%
% @param: DATA : case number or name for sequence to run
%----------------------------------------------------------------------------%
function [out, etc] = cdov_first_last(DATA, model, versiontype)
%----------------------------------------------------------------------------%
% some startup + parameters
%----------------------------------------------------------------------------%
setup;

%----------------------------------------------------------------------------%
% Settings
%----------------------------------------------------------------------------%
if ~exist('versiontype', 'var');
  % versiontype = 'do_frame_0';
  versiontype = 'dov_frame_ayvaciuv_70_fxf_pfg0p001';
end
BEGIN = 1;
%----------------------------------------------------------------------------%
% DATA PATH
%----------------------------------------------------------------------------%
if isnumeric(DATA); % btay datapath
  paramsIn = struct();
  paramsIn.expName = 'combine';
  [p, params] = ttaobdData(DATA, paramsIn);
  p.expName = paramsIn.expName;
  img_path = p.dpath;
  flow_path = p.uvpath;
  seq = p.seq;

  % how long to go for
  imgNameSearchStr = regexprep(p.imgNameStr, '\%\d*d', '*');
  files = dir(fullfile(img_path, sprintf(imgNameSearchStr, seq)));
  if isempty(files); files = dir([img_path '/*.' p.ext]); end
  uvNameSearchStr  = regexprep(p.uvNameStr,  '\%\d*d', '*');
  flow_files = dir(fullfile(flow_path, sprintf(uvNameSearchStr, seq)));
  if isempty(flow_files); flow_files = dir([flow_path '/*.mat']); end
  
else % vasiliy datapath
  [seq, flow_path, img_path, ~, extension] = dataPaths(DATA);
  files = dir([img_path '/*.' extension]);
  flow_files = dir([flow_path '/*.mat']);
end

% checks for rsz files
[files, ~] = purge_files(files, 'rsz');
T = length(files);

%----------------------------------------------------------------------------%
% algorithm settings
%----------------------------------------------------------------------------%
params = tao_params_default();

%==============================================================================
% use this to overwrite any variables used in here before running
% NOTE: all variables should be set above this point
if exist('PKG','var') && isstruct(PKG); params = structmerge(params, PKG); end;
%==============================================================================
%--------------------------------------------------------------------------
% output files
%--------------------------------------------------------------------------
if exist('params', 'var') && isfield(params, 'outpath');
  outpath = params.outpath;
else
  if ~exist('model', 'var'); model = 'full'; end;
  if strcmp(model, 'fxf');
    outpath = ['/plot/btay/projects/detachable/cdov/week17-cdov-r0.1-fxf'];
  else
    outpath = ['/plot/btay/projects/detachable/cdov/week17-cdov-r0.1'];
  end
end

% add sequence name to outpath
outpath = [outpath, '/', seq];
if ~exist(outpath,'dir'), mkdir(outpath); end;
if ~exist('nameStr', 'var'); nameStr = '%s_%06d'; end;

%==============================================================================
% run the sequence
%==============================================================================
BEGIN = BEGIN + 1;
FINISH = T - 1;
for k = [BEGIN, FINISH];
  fprintf('=========================== time: %d =============\n', k);
  %------------------------------------------------------------------
  % image data
  %------------------------------------------------------------------
  I1 = imread(fullfile(img_path, files(k).name));
  I1_bflt = recursive_bf_mex(I1, 0.01, 0.1, 1, 5);
  i1_bflt = im2double(I1_bflt);

  %------------------------------------------------------------------
  % flow + occlusions
  %------------------------------------------------------------------
  if (k == 2);
    uvb_rev_file = fullfile(flow_path, flow_files(k - 1).name);
    uvb_rev_data = load(uvb_rev_file);
    uvb_rev = uvb_rev_data.uvf;
    occb_rev = uvb_rev_data.ef;
    
    occb_rev(isnan(occb_rev) | isinf(occb_rev)) = 0.0;
    occb_rev_prob = occ_res_to_prob(occb_rev, true);
    uvb_rev_bflt = recursive_bf_mex(uvb_rev, 0.05, 0.004, 0, 10);
  elseif (k == FINISH);
    uvf_rev_file = fullfile(flow_path, flow_files(k + 1).name);
    uvf_rev_data = load(uvf_rev_file);
    uvf_rev = uvf_rev_data.uvb;
    occf_rev = uvf_rev_data.eb;

    occf_rev(isnan(occf_rev) | isinf(occf_rev)) = 0.0;
    occf_rev_prob = occ_res_to_prob(occf_rev, true);
    uvf_rev_bflt = recursive_bf_mex(uvf_rev, 0.05, 0.004, 0, 10);
  end
  
  %------------------------------------------------------------------
  % load old layers
  %
  % TODO: should use version type to not grab the wrong ones
  %------------------------------------------------------------------
  outResFile = fullfile(outpath, sprintf(nameStr, seq, k));
  lay_files = dir([outResFile, '*lay*.mat']);
  data = load(fullfile(outpath, lay_files(1).name), 'lay');
  lay_t0 = data.lay;
  
  %------------------------------------------------------------------
  % 2nd frame and last frame to label
  %------------------------------------------------------------------
  if (k == 2);
    outRes1File = fullfile(outpath, sprintf(nameStr, seq, k-1));
    % lay = warp_output_layers(lay_t0, uvb_rev_bflt, occb_rev_prob > OCCPROB);
    lay = warp_output_layers(lay_t0, uvb_rev_bflt, ...
      occb_rev_prob > params.OCCPROBLAYERTHRESH);
    
  elseif (k == FINISH);
    outRes1File = fullfile(outpath, sprintf(nameStr, seq, k+1));
    lay = warp_output_layers(lay_t0, uvf_rev_bflt, ...
      occf_rev_prob > params.OCCPROBLAYERTHRESH);
  end
  lay = postfilter_layers(lay, i1_bflt, 0.5);
  save([outRes1File, '_lay_', versiontype, '.mat'], 'lay');
end
unix(sprintf('touch %s/%s.done', outpath, sprintf(nameStr, seq, k+1)));
exit;
end
