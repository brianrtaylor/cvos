%--------------------------------------------------------------------
% make_unary_constraints
%
% makes unary constraint map to be used by optimization
%--------------------------------------------------------------------
% function [map] = make_unary_constraints(past, ...
%   constraint_uno_img_now_b, constraint_uno_img_now_f, ...
%   constraint_weights_uno_img_now_b, constraint_weights_uno_img_now_f)
% function [map] = make_unary_constraints(past, img_b, img_f, wx_l, wy_l, opts)
function [map, past_map] = make_unary_constraints(past, img_b, img_f, opts)

v2struct(opts);
% if ~exist('wx_l', 'var'); wx_l = 0; end;
% if ~exist('wy_l', 'var'); wy_l = 0; end;
[rows, cols, ~] = size(img_b);

% TODO: add the following
% * removal of unary constraints that don't satisfy the final result, 
%   so another \xi persay (\xi_unary)
% * UNARY_CONSTRAINT_MEMORY_FACTOR
% * DO_UNARY_CONSTRAINT_WARP_WEIGHT : thus if the frame / etc is moving,
%     eventually we remove the 
% * don't apply the diminishing to those cues that satisfy the layer results

% past.layers to validate / diminish the past.unary_constraints_map

% mask_past_layer = reshape(min(0.0, past.layers), [rows, cols]) - wx_l - wy_l;
mask_past_layer = zeros(rows, cols);
if ~isempty(past.layers);
  mask_past_layer = reshape(min(1.0, past.layers), [rows, cols]);
end
blob = strel('disk', 5);
smallblob = strel('disk', 1);

past_map = past.unary_constraints_map .* max(mask_past_layer, ...
  (1 - UNARY_CONSTRAINT_MEMORY_FACTOR));

if DO_UNARY_CONSTRAINT_WARP_WEIGHT && ~isempty(past.uvf_rev);
  uv_mag = sqrt(sum(past.uvf_rev .^ 2, 3));
  uv_mag_avg = mean(uv_mag(uv_mag > 1.0));
  
  w_uv_mag = exp( -uv_mag ./ max(WARP_SAFESPEEDSQUARED, 2 * uv_mag_avg) );
  % weighting for inaccurate motion
  past_map = past_map .* max(mask_past_layer, w_uv_mag);
  % past_map = past_map .* w_uv_mag;
end

past_map = max(0.0, past_map - imdilate(img_f, smallblob));

% forward adds, backward removes
map = max(0.0, past_map + img_f - imdilate(img_b, blob));

if VIS < 100;
  maxval = 5.0;
  fig(415); imagesc([past_map, map, mask_past_layer * maxval], ...
    [0, maxval]); colorbar; title('occlusion_unary_constraint');
  drawnow;
end
end
