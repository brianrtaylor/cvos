%------------------------------------------------------------------
% tao weights: computes weights step for framework
% 
% @param: opts : containts : edge_model 
%------------------------------------------------------------------
% function [weights, wx_l, wy_l] = propagate_tao_weights_old_frame( ...
%   weights_now, I1, uvb, uvb_cbf, uvf, uvf_cbf, uvf_old, ...
%   constraints_now_b, constraint_weights_now_b, ...
%   constraints_now_f, constraint_weights_now_f, ...
%   layers, prev_layers, past_weights, occb_cbf_prob, ...
%   Dx, Dy, dx_inds, dy_inds, past, k, opts)
function [unity, wx_l, wy_l, past] = tao_unity_prior( ...
  past, uvb, uvb_cbf, wx_l, wy_l, occb_mask, k, opts)

%------------------------------------------------------------------
% preliminaries
%------------------------------------------------------------------
% get options
ERODE_SZ = 3.0;
v2struct(opts);

[rows, cols, ~] = size(occb_mask);
imsize = [rows, cols];
n = rows * cols;

g_tiny = fspecial('gauss', 5, 1);

mag_uvb_squared = sum(uvb_cbf .^ 2, 3);
mag_uvb = sqrt(mag_uvb_squared);

%------------------------------------------------------------------
% weights from prior frame
%------------------------------------------------------------------
wx = past.unity(:, :, 1); % reshape(past.unity(1:n), imsize);
wy = past.unity(:, :, 2); % reshape(past.unity((n+1):(2*n)), imsize);

% using layers
unity_mask = max(0.0, min(past.layers, 1.0) - wx_l - wy_l);
wx = wx .* (1.0 - UNITY_MEMORY_FACTOR) + unity_mask;
wy = wy .* (1.0 - UNITY_MEMORY_FACTOR) + unity_mask;

% diminish weights due to warp: faster motion means less trust in weights
if DO_UNITY_WARP_WEIGHT; % increase in wx_time = trust weight less (por)
  % TODO: if I uncomment these two lines below to make this more "proper"
  % in my opinion, then I lose the rear car window in cars1. hmmm :/. so
  % far I leave it commented out to keep that window, and keep that unity
  % prior high
  %   wx = wx .* past.w_warp_uvf; % used to be commented out 20141027
  %   wy = wy .* past.w_warp_uvf; % used to be commented out 20141027
  wx = wx .* past.w_fg_uvf;
  wy = wy .* past.w_fg_uvf;
end

% unity mask t0
wx_l_t0 = utils_warp_image(wx_l, uvb);
wy_l_t0 = utils_warp_image(wy_l, uvb);
unity_mask_t0 = max(0.0, min(past.t0.layers, 1.0) - wx_l_t0 - wy_l_t0);
occb_mask = max(0.0, occb_mask - unity_mask_t0 * ...
  OCC_INSIDE_UNITY_MEMORY_FACTOR);

% create the unity prior layers things
wx_t0 = zeros(imsize);
wy_t0 = zeros(imsize);
n_past_layers = max(past.layers(:));

% does this do anything for kk = 0? <-- check it
for kk = 0:n_past_layers;
  layer_indicator = abs(past.layers - kk) < 0.1;
  wx_ul = wx .* layer_indicator;
  wy_ul = wy .* layer_indicator;

  % start old normal stuff
  wx_ul_t0 = utils_warp_image(wx_ul, uvb);
  wy_ul_t0 = utils_warp_image(wy_ul, uvb);

  wx_ul_t0(isnan(wx_ul_t0) | isinf(wx_ul_t0)) = 0;
  wy_ul_t0(isnan(wy_ul_t0) | isinf(wy_ul_t0)) = 0;

  % erode edges a bit based on flow magnitude in that region
  [wx_ul_t0, wx_msfm] = remove_occb_mask(wx_ul_t0, occb_mask, ERODE_SZ, 1);
  wx_diminish = 1 - exp(-max(0.0, wx_msfm - 2*ERODE_SZ) ./ max(1.0, mag_uvb));
  wx_ul_t0 = wx_ul_t0 .* wx_diminish;
  [wy_ul_t0, wy_msfm] = remove_occb_mask(wy_ul_t0, occb_mask, ERODE_SZ, 1);
  wy_diminish = 1 - exp(-max(0.0, wy_msfm - 2*ERODE_SZ) ./ max(1.0, mag_uvb));
  wy_ul_t0 = wy_ul_t0 .* wy_diminish;
  
  wx_ul_t0 = imfilter(wx_ul_t0, g_tiny);
  wy_ul_t0 = imfilter(wy_ul_t0, g_tiny);

  wx_t0 = wx_t0 + wx_ul_t0;
  wy_t0 = wy_t0 + wy_ul_t0;
end

wx_t0(isnan(wx_t0) | isinf(wx_t0)) = 0;
wy_t0(isnan(wy_t0) | isinf(wy_t0)) = 0;
past.t0.unity = cat(3, wx_t0, wy_t0); % not great

unity = cat(3, wx_t0, wy_t0);

%------------------------------------------------------------------
% weights from prior frame visualization
%------------------------------------------------------------------
% diagnostics
if VIS < 160;
wait_img = [[wx, wy, min(wx, wy)]; ...
  [past.w_warp_uvf, past.w_fg_uvf, past.w_warp_uvf .* past.w_fg_uvf]; ...
  [wx_t0, wy_t0, max(wx_t0, wy_t0)]];
fig(4); clf; imagesc(wait_img, [0, max(1, max(wait_img(:)))]);
end

if VIS < 180;
  wwp = vis_weights_occlusions(past.t0.unity, [], imsize, []);
  ww = vis_weights_occlusions(unity, [], imsize, []);
  
  fig(119); clf;
  vl_tightsubplot(2,2,1);
  imagesc(wwp); axt; axis off; title('old unity prior (b4 warp)');
  vl_tightsubplot(2,2,2);
  imagesc(ww); axt; axis off; title('unity prior');
  vl_tightsubplot(2,2,3);
  imagesc(max(wx, wy)); axt; axis off; title(sprintf('wxy: %d', k));
  vl_tightsubplot(2,2,4);
  imagesc(max(wx_t0, wy_t0)); axt; axis off; title('wxy warped');
  drawnow;
end
end
