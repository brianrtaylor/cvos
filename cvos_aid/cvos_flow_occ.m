%-----------------------------------------------------------------------------
% cvos_flow_occ
%
% function to handle all the loading of optical flow from file and performing
% cross bilateral filtering for extrapolating flow in the occluded region
%
% @return: uv{b,f}: optical flow from t to t-1 and t+1
% @return: uv{b,f}_cbf: uv{b,f} after cross-bilateral filter
% @return: uv{b,f}_rev: optical flow from t-1 and t+1 to t
% @return: occb{,_cbf,_cbf_prob}: occluded region going from t to t-1, post-
%   cross bilateral filtering, and after being converted to a probability
% @return: occf{,_cbf,_cbf_prob}: occluded region going from t to t+1, post-
%   cross bilateral filtering, and after being converted to a probability
% @return: occ{b,f}_rev: occluded region going from t-1 and t+1 to t
% @return: occ{b,f}_rev_prob: occ{b,f}_rev converted to a probability
%
% @param: uvpath: path to optical flow files
% @param: flow_files: struct with flow file names and more info
% @param: seq: sequence name of interest to process
% @param: k: current frame number
% @param: T: maximum frames to process
% @param: i0, i1, i2: images t-1, t, and t+1 with range [0,1]
% @param: I0, I1, I2: images t-1, t, and t+1 with range [0,255]
% @param: past: struct with past information so we can reuse old loaded flows
% @param: opts: struct with parameters
%-----------------------------------------------------------------------------
function [uvb, uvb_cbf, uvf, uvf_cbf, uvb_rev, uvf_rev, ...
  occb, occb_cbf, occb_cbf_prob, occf, occf_cbf, occf_cbf_prob, ...
  occb_rev, occf_rev, occb_rev_prob, occf_rev_prob] = cvos_flow_occ( ...
  uvpath, flow_files, seq, k, T, i0, i1, i2, I0, I1, I2, past, opts)

DO_CROSSBILATERALFILTERFLOW = false;
v2struct(opts);
num_files = length(flow_files);

%-----------------------------------------------------------------------------
% loading flow files
%-----------------------------------------------------------------------------
tic;
if CAUSAL && DO_FORBACKCAUSAL && k >= num_files;
  kk = 2 * num_files - k;

  uv_file      = fullfile(uvpath, flow_files(kk    ).name); 
  uvf_rev_file = fullfile(uvpath, flow_files(kk - 1).name); 
  
  uv_data = load(uv_file);
  uvf_rev_data = load(uvf_rev_file);
  
  uvb_rev = past.uvf;
  uvb     = past.uvf_rev;
  uvf     = uv_data.uvb;
  uvf_rev = uvf_rev_data.uvf;
  
  occb_rev  = past.occf_cbf;
  occb      = past.occf_rev;
  occf      = uv_data.eb;
  occf_rev  = uvf_rev_data.ef;
  
elseif (k == 1) || (k == T && opts.DO_FORBACKCAUSAL); % first frame
  DO_CROSSBILATERALFILTERFLOW = false;
  uv_file      = fullfile(uvpath, flow_files(k    ).name); 
  uvf_rev_file = fullfile(uvpath, flow_files(k + 1).name); 
  
  uv_data = load(uv_file);
  uvf_rev_data = load(uvf_rev_file);

  uvf       = uv_data.uvf;
  uvf_rev   = uvf_rev_data.uvb;
  occf      = uv_data.ef;
  occf_rev  = uvf_rev_data.eb;
  
  [rows, cols, ~] = size(uvf);
  uvsize = [rows, cols, 2];
  imsize = [rows, cols];
  
  uvb_rev   = zeros(uvsize);
  uvb       = zeros(uvsize);
  occb_rev  = zeros(imsize);
  occb      = zeros(imsize);    
  
elseif k == num_files && ~opts.DO_FORBACKCAUSAL; % last frame
  DO_CROSSBILATERALFILTERFLOW = false;
  uvb_rev_file = fullfile(uvpath, flow_files(k - 1).name); 
  uv_file      = fullfile(uvpath, flow_files(k    ).name); 
  
  uv_data = load(uv_file);
  uvb_rev_data = load(uvb_rev_file);

  uvb       = uv_data.uvb;
  uvb_rev   = uvb_rev_data.uvf;
  occb_rev  = uvb_rev_data.ef;
  occb      = uv_data.eb;    

  [rows, cols, ~] = size(uvb);
  uvsize = [rows, cols, 2];
  imsize = [rows, cols];    

  uvf       = zeros(uvsize);
  uvf_rev   = zeros(uvsize);
  occf      = zeros(imsize);
  occf_rev  = zeros(imsize);

elseif CAUSAL && ~isempty(past.uvf); % causal and not first frame
  uv_file      = fullfile(uvpath, flow_files(k    ).name); 
  uvf_rev_file = fullfile(uvpath, flow_files(k + 1).name); 
  
  uv_data = load(uv_file);
  uvf_rev_data = load(uvf_rev_file);
  
  uvb_rev = past.uvf;
  uvb     = past.uvf_rev;
  uvf     = uv_data.uvf;
  uvf_rev = uvf_rev_data.uvb;
  
  occb_rev  = past.occf_cbf;
  occb      = past.occf_rev;
  occf      = uv_data.ef;
  occf_rev  = uvf_rev_data.eb;
  
else % frame 2 (first useable frame)
  uvb_rev_file = fullfile(uvpath, flow_files(k - 1).name); 
  uv_file      = fullfile(uvpath, flow_files(k    ).name); 
  uvf_rev_file = fullfile(uvpath, flow_files(k + 1).name); 
  
  uv_data = load(uv_file);
  uvf_rev_data = load(uvf_rev_file);
  uvb_rev_data = load(uvb_rev_file);
  
  uvb_rev = uvb_rev_data.uvf;
  uvb     = uv_data.uvb;
  uvf     = uv_data.uvf;
  uvf_rev = uvf_rev_data.uvb;
  
  occb_rev  = uvb_rev_data.ef;
  occb      = uv_data.eb;
  occf      = uv_data.ef;
  occf_rev  = uvf_rev_data.eb;
end
toc; % time it takes to load data

occb_rev(isnan(occb_rev) | isinf(occb_rev)) = 0.0;
occb    (isnan(occb    ) | isinf(occb    )) = 0.0;
occf    (isnan(occf    ) | isinf(occf    )) = 0.0;
occf_rev(isnan(occf_rev) | isinf(occf_rev)) = 0.0;

occb_rev_prob = occ_res_to_prob(occb_rev, true);
occb_prob     = occ_res_to_prob(occb    , true);
occf_prob     = occ_res_to_prob(occf    , true);
occf_rev_prob = occ_res_to_prob(occf_rev, true);

%-----------------------------------------------------------------------------
% flow
%-----------------------------------------------------------------------------
if DO_CROSSBILATERALFILTERFLOW;
  [uvb_cbf_, uvf_cbf_, occb_cbf_, occf_cbf_, occb_cbf_prob_, ...
    occf_cbf_prob_] = cross_bilateral_filter_step_test_fullimg( ...
    double(occb), double(occf), double(uvb), double(uvf), ...
    double(uvb_rev), double(uvf_rev), I0, I1, I2, m, VIS);
  [uvb_cbf, uvf_cbf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob] ...
    = cross_bilateral_filter_step_test_occ(double(occb_cbf_prob_), ...
    double(occf_cbf_prob_), double(uvb), double(uvb_cbf_), double(uvf), ...
    double(uvf_cbf_), double(occb_rev), double(occf_rev), double(uvb_rev), ...
    double(uvf_rev), I0, I1, I2, m, VIS);
else
  uvb_cbf = uvb;
  uvf_cbf = uvf;
  occb_cbf = occb;
  occf_cbf = occf;
  occb_cbf_prob = occb_prob;
  occf_cbf_prob = occf_prob;
end

if VIS < 150;
  fig(111); imagex([[uvb, uvf]; [uvb_cbf, uvf_cbf]; [uvb_rev, uvf_rev]]);
  title('cvos flow (cbf)'); colorbar; drawnow;
  fig(112); imagesc([[occb_prob, occf_prob]; [occb_cbf_prob, ...
    occf_cbf_prob]; [occb_cbf_prob, occf_cbf_prob] > OCCPROB]);
  title('cvos occ (cbf)'); colorbar; drawnow;
end
end
