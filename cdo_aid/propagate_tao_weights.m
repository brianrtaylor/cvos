%------------------------------------------------------------------
% tao weights: computes weights step for framework
% 
% @param: opts : contraints : edge_model 
%------------------------------------------------------------------
% function [weights, wx_l, wy_l, past] = propagate_tao_weights( ...
%   weights_now, uvb, past, occb_mask, ...
%   Dx, Dy, k, opts, object_mean_uv_map)
%------------------------------------------------------------------
function [weights, wx_l, wy_l, past] = propagate_tao_weights( ...
  weights_now, uvb, past, occb_mask, Dx, Dy, k, opts)

%------------------------------------------------------------------
% preliminaries
%------------------------------------------------------------------
v2struct(opts);
[rows, cols, ~] = size(uvb);
imsize = [rows, cols];
n = rows * cols;

%------------------------------------------------------------------
% weights from prior frame
%------------------------------------------------------------------
% old weights (positive)
wx = past.weights(:, :, 1);
wy = past.weights(:, :, 2);

% using layers
dx_l = Dx * past.layers(:);
dy_l = Dy * past.layers(:);
wx_l = reshape(dx_l, imsize);
wy_l = reshape(dy_l, imsize);
wx_l = min(abs(wx_l), 1.0); % 1 - change in layer (used), 0 - same layer
wy_l = min(abs(wy_l), 1.0);

% diminish weights due to age
wx = wx * (1 - WEIGHT_MEMORY_FACTOR);
wy = wy * (1 - WEIGHT_MEMORY_FACTOR);

% diminish weights due to warp: faster motion means less trust in weights
if DO_WEIGHT_WARP_WEIGHT; % increase in wx_time = trust weight less (por)
  wx = wx .* past.w_warp_uvf;
  wy = wy .* past.w_warp_uvf;
  wx = wx .* past.w_fg_uvf;
  wy = wy .* past.w_fg_uvf;
end
  
% prior weights warped in for visualization later
% TODO: can combine these two warpings into 1 by adding wx and wx_l first
wx_t0 = utils_warp_image(wx, uvb);
wy_t0 = utils_warp_image(wy, uvb);

wx_t0 = remove_occb_mask(wx_t0, occb_mask, 0, 1);
wy_t0 = remove_occb_mask(wy_t0, occb_mask, 0, 1);

wx_t0_l = utils_warp_image(wx_l, uvb);
wy_t0_l = utils_warp_image(wy_l, uvb);

%------------------------------------------------------------------
% if edge was used, boost it with w{x,y}_t0_l
%------------------------------------------------------------------
wx_t0_f = wx_t0 + wx_t0_l;
wy_t0_f = wy_t0 + wy_t0_l;
wx_t0_f(isnan(wx_t0_f) | isinf(wx_t0_f)) = 0;
wy_t0_f(isnan(wy_t0_f) | isinf(wy_t0_f)) = 0;

past.t0.weights = cat(3, wx_t0_f, wy_t0_f);


%------------------------------------------------------------------
% weights edge weight
%------------------------------------------------------------------
wnx_img = weights_now(:, :, 1);
wny_img = weights_now(:, :, 2);
  %   % spreading / contracting
% %   % wnx_img = imfilter(wnx_img, g);
% %   % wny_img = imfilter(wny_img, g);
% %   weights_now = cat(1, wnx_img(:), wny_img(:));

% assert(all(vec(wx_t0_f > 0)), 'wx_t0_f < 0');
% assert(all(vec(wy_t0_f > 0)), 'wx_t0_f < 0');

wfx_img = wnx_img + wx_t0_f;
wfy_img = wny_img + wy_t0_f;
weights = cat(3, wfx_img, wfy_img);

%------------------------------------------------------------------
% visualization
%------------------------------------------------------------------
if VIS < 180;
  % diagnostics
  % wait_img = [[wx, wy, max(wx, wy)]; ...
  %   [wx_l   , wy_l   , max(wx_l   , wy_l   )]; ...
  %   [wx_l+wx, wy_l+wy, max(wx_l+wx, wy_l+wy)]; ...
  %   [wnx_img, wny_img, max(wnx_img, wny_img)]; ...
  %   [wfx_img, wfy_img, max(wfx_img, wfy_img)]];
  wait_img = [[wx_t0, wy_t0, max(wx_t0, wy_t0)]; ...
    [wx_t0_l, wy_t0_l, max(wx_t0_l, wy_t0_l)]; ...
    [wx_t0_f, wy_t0_f, max(wx_t0_f, wy_t0_f)]; ...
    [wnx_img, wny_img, max(wnx_img, wny_img)]; ...
    [wfx_img, wfy_img, max(wfx_img, wfy_img)]];
  fig(3); clf; sc(wait_img, 'jet', [0, max(1, max(wait_img(:)))]);
end

if VIS < 250;
  ww  = vis_weights_occlusions(max(0.0, 1.0 - weights_now), [], imsize, []);
  ww_ = vis_weights_occlusions(max(0.0, 1.0 - weights    ), [], imsize, []);
  wwp  = vis_weights_occlusions(max(0.0, 1.0 ...
    - cat(1, wx(:)   , wy(:)))   , [], imsize, []);
  wwpw = vis_weights_occlusions(max(0.0, 1.0 ...
    - cat(1, wx_t0(:), wx_t0(:))), [], imsize, []);

  fig(118); clf;
  vl_tightsubplot(3,2,1);
  % imagesc( min(wx, wy) ); axt; title('old weights b4 warp');
  imagesc(wwp); axt; axis off; title('old weights b4 warp');
  vl_tightsubplot(3,2,2);
  % imagesc( min(wx_warped, wy_warped) ); axt; title('old weights warped');
  imagesc(wwpw); axt; axis off; title('old weights warped');
  
  vl_tightsubplot(3,2,3);
  imagesc(min(1.0-wx, 1.0-wy)); axt; axis off;
  title(sprintf('wxy: %d', k));
  
  vl_tightsubplot(3,2,4);
  imagesc(min(1.0-wx_t0, 1.0-wy_t0)); axt; axis off;
  title('wxy warped');
  
  vl_tightsubplot(3,2,5);
  imagesc(ww); axt; axis off; title('current weights');
  vl_tightsubplot(3,2,6);
  imagesc(ww_); axt; axis off; title('final weights');
  drawnow;
end

%------------------------------------------------------------------
% quick error check
%------------------------------------------------------------------
if any(vec(weights(:) < weights_now(:)));
  fprintf('%s: weights got weaker?!\n', mfilename);
end
end
