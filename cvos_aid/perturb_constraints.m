%-----------------------------------------------------------------------------
% perturb_constraints
%
% obtain shape model for local shape classifier from segmentation from the 
% previous frame warped into the current frame
%
% @return: constraints_out: resulting perturbed constraints
% @return: gmm_weights: weights computed from the gmm used for perturbation
%
% @param: constraints: input constraints naively warped by optical flow
% @param: constraint_groups: constraints are perturbed within local groups
%   identified by these group ids
% @param: local_gmm_{f,b}g: gmm parameters 
% @param: I0: input rgb image
% @param: weights: segmentation weights computed on current frame
% @param: params: parameters struct
%-----------------------------------------------------------------------------
function [constraints_out, gmm_weights] = perturb_constraints(constraints, ...
  constraint_groups, local_gmm_fg, local_gmm_bg, I0, weights, params)

if ~isa(I0,'double'); I0 = im2double(I0); end;

imsize = [size(I0,1), size(I0,2)];
W = 0.5 * sum(weights, 3);
W_binary = W <= 0.75;
SpeedImage = ones(imsize);
[ii,jj] = find(W_binary);
dist_edge = msfm(SpeedImage, [ii,jj]', true, true);
d_edge = exp(-0.5 * dist_edge);
clear SpeedImage W W_binary

[YY, XX] = ind2sub(imsize, constraints(:,1) );
pts_occr = [XX(:), YY(:)];   
[YY, XX] = ind2sub(imsize, constraints(:,2) );
pts_occd = [XX(:), YY(:)];

pts_max_separation = params.MAXCONSTRAINTDIST * ones(size(pts_occr, 1), 1);

[pts_occd_out, pts_occr_out, gmm_weights] = perturb_constraints_mex( ...
  I0, d_edge, pts_occd, pts_occr, pts_max_separation, constraint_groups, ...
  local_gmm_fg, local_gmm_bg, params);
try
  idx_occd = sub2ind(imsize, pts_occd_out(:,2), pts_occd_out(:,1));
  idx_occr = sub2ind(imsize, pts_occr_out(:,2), pts_occr_out(:,1));
catch
  disp('failed'); 
end

gmm_weights(isnan(gmm_weights) | isinf(gmm_weights)) = 0.0;
constraints_out = [idx_occr(:), idx_occd(:)];
end
