%-----------------------------------------------------------------------------
% vis_weights_occlusions
%
% visualizes the occlusions on top of the weights image
%
% @return: W: output image with cues 
% @param: weights (MxNx2): weights image to draw occlusion cues on top of
% @param: constraints: occluder/occluded pairs to draw on image
% @param: imsize: image dimensions (rows, cols)
% @param: constraint_weights: associated weights for each of the constraints
%-----------------------------------------------------------------------------
function W = vis_weights_occlusions(weights, constraints, imsize, ...
  constraint_weights)
nconstraints = size(constraints, 1);
if ~exist('constraint_weights', 'var'); 
  constraint_weights = ones(nconstraints, 1);
end

constraint_weights = min(1.0, constraint_weights);
weights = min(1.0, weights);

nnodes = numel(weights) / 2;
W = reshape(min(weights(1:nnodes), weights(nnodes + 1:end)), imsize);
W = repmat(W, [1 1 3]);
W = max(min(W, 1.0), 0.0);

z0 = zeros(imsize);
if ~isempty(constraints);
  % to ensure that the largest value is placed in a pixel (assuming 
  % last is set in underlying matlab functions)
  [constraint_weights_, ind] = sort(constraint_weights);
  constraints_ = constraints(ind, :);
  
  occd = z0;
  occd(constraints_(:, 2))=constraint_weights_;
  occd = cat(3, occd, z0, z0); % red
  
  occr = z0;
  occr(constraints_(:,1))=constraint_weights_;
  occr = cat(3, occr, occr, z0); % yellow

  W = 0.3 * W + 0.7 * (W .* repmat(occd(:, :, 1) == 0, [1 1 3]) + occd);
  W = 0.3 * W + 0.7 * (W .* repmat(occr(:, :, 1) == 0, [1 1 3]) + occr);
end
end
