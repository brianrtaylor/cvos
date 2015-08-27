%-----------------------------------------------------------------------------
% make_cvos_weights_current_frame
%
% computes pixelwise weights for layer segmentation on the current frame
%
% @param: Ilab: image to compute weights from
% @param: uv{b,f}: backward and forward optical flow from t to t-1 and t+1
% @param: d{x,y}_inds: indices for edges in x and y directions
% @param: E: edge strength
% @param: opts: structure of parameters
%-----------------------------------------------------------------------------
function [weights, weights_inds] = make_cvos_weights_current_frame( ...
  Ilab, uvb, uvf, dx_inds, dy_inds, E, opts)
v2struct(opts);

weights = make_pixelwise_weights([dx_inds; dy_inds], Ilab, uvb, uvf, E, ...
  relative_weights, opts);

weights_sigmoid = sigmoid(weights, SIGMOID_MEAN, SIGMOID_SCALE);
weights = weights_sigmoid;
weights( weights < WEIGHTS_LOW_CUTOFF ) = 0;

weights = 1 - weights;
weights_inds = [dx_inds; dy_inds];
end
