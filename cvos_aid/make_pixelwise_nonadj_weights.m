%-----------------------------------------------------------------------------
% make_pixelwise_nonadj_weights
%
% computes pixelwise weights on the current frame using intensity, flow, and
% edge information
%
% @return: w: final computed weights
% @return: w_{div,img,uvb,uvf,edgecross}: the constraint weights computed 
%   from various input feature channels (the last of which depends on whether
%   or not the constraint straddles a strong edge response)
% @param: idx: occluder-occluded constraints
% @param: I: current image frame
% @param: uv{b,f}_ (MxNx2): backward and forward flow from t to t-1 and t+1
% @param: gpb: edge map strength
% @param: weights: initial weights computed elsewhere, which we multiply by w
% @param: DIVFLAG: whether to use divergence to affect constraint weight
%-----------------------------------------------------------------------------
function [w, w_div, w_img, w_uvb, w_uvf, w_edgecross] ...
  = make_pixelwise_nonadj_weights(idx, I, uvb_, uvf_, gpb, weights, DIVFLAG)
if ~exist('DIVFLAG', 'var'); DIVFLAG = false; end;
weights = weights / sum(weights);
[rows, cols, chan] = size(I);
I = reshape(I, [numel(I) / 3, 3]);
uvf = reshape(uvf_, [numel(uvf_) / 2, 2]);
uvb = reshape(uvb_, [numel(uvb_) / 2, 2]);
idx1 = idx(:, 1);
idx2 = idx(:, 2);

% intensity weight
if weights(1) > 0;
  w_img = sum((I(idx1, :) - I(idx2, :)) .^ 2, 2);
  w_img = exp(-w_img / mean(w_img));
else
  w_img = zeros(length(idx1), 1);
end

if weights(2) > 0;
  % flow weight (forward)
  w_uvf = sum((uvf(idx1, :) - uvf(idx2, :)) .^ 2, 2);
  w_uvf = exp(-w_uvf / mean(w_uvf));
  % flow weight (backward)
  w_uvb = sum((uvb(idx1, :) - uvb(idx2, :)) .^ 2, 2);
  w_uvb = exp(-w_uvb / mean(w_uvb));
else
  w_uvb = zeros(length(idx1), 1);
  w_uvf = zeros(length(idx1), 1);        
end

% gpb weight
if weights(3) > 0;
  gpb = gpb(:);
  w_gpb = min(gpb(idx1), gpb(idx2));
else
  w_gpb = zeros(length(idx1), 1);
end

% edge crossing weight
if weights(4) > 0;
  w_edgecross = occweight_edge_crossing(double(gpb), [idx1, idx2]);
else
  w_edgecross = zeros( length(idx1), 1);
end

w = w_img * weights(1) + (w_uvf + w_uvb) * weights(2)/2 ...
  + w_gpb * weights(3) + w_edgecross * weights(4);

w(isnan(w) | isinf(w)) = 0.0;

%-----------------------------------------------------------------------------
% divergence scaling factor
%-----------------------------------------------------------------------------
w_div = [];
if DIVFLAG;
  [y1, x1] = ind2sub([rows, cols], idx1);
  [y2, x2] = ind2sub([rows, cols], idx2);
  
  vector_flow_uvf = uvf(idx1,:) - uvf(idx2,:);
  vector_flow_uvb = uvb(idx1,:) - uvb(idx2,:);
  vector_loc = ([x1, y1] - [x2, y2]);
  vector_loc_denom = repmat(sqrt(sum(vector_loc .^ 2, 2)), [1, 2]);
  vector_loc = vector_loc ./ vector_loc_denom;
  
  divf = sum(vector_flow_uvf .* vector_loc, 2);
  divb = sum(vector_flow_uvb .* vector_loc, 2);

  DIVDENOM = 0.5;
  w_div = max(0, 1 - exp((divf - divb) / DIVDENOM));
  w_div(isnan(w_div) | isinf(w_div)) = 0.0;
end
end
