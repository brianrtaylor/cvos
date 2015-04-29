%------------------------------------------------------------------
% tao weights: computes weights step for framework
%
% @param: I1 : image to compute weights from
% @param: ...
% @param: opts : containts following
%   * edge_model
%   * SIGMOID_MEAN
%   * SIGMOID_SCALE
%   * relative_weights
%   * relative_weights_constraints
%------------------------------------------------------------------
function [weights, weights_inds] = make_tao_weights_current_frame( ...
  Ilab, uvb, uvf, dx_inds, dy_inds, E, opts)
v2struct(opts);

%------------------------------------------------------------------
% weights current frame
%------------------------------------------------------------------ 
weights = make_pixelwise_weights2([dx_inds; dy_inds], ...
  Ilab, uvb, uvf, E, relative_weights, opts);

weights_sigmoid = sigmoid(weights, SIGMOID_MEAN, SIGMOID_SCALE);
weights = weights_sigmoid;
weights( weights < WEIGHTS_LOW_CUTOFF ) = 0;

% hack to make higher values wanted
weights = 1 - weights;
% weights = min(weights_sigmoid, weights); % vasiliy hack
weights_inds = [dx_inds; dy_inds];
end
