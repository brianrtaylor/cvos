%-----------------------------------------------------------------------------
% make_cvos_constraints
%
% constructs a list of occluder-occluded constraints and weights from the 
% set of occluded pixels and the forward and reverse flow to another frame
%
% @param: occ_prob: probability of occlusion at each pixel
% @param: uv: forward warping from current frame to frame X
% @param: uv_rev: reverse warping from frame X to current frame
% @param: OCCPROB: threshold for occ_prob being occluded
% @param: MINCONSTRAINTDIST: minimum separation between occluder and occluded
%   pixel for constraint to be valid
% @param: VIS: visualization level parameter
% @param: fig_idx: figure number for visualization
%-----------------------------------------------------------------------------
function [constraints, constraint_weights, bad] = make_cvos_constraints( ...
  occ_prob, uv, uv_rev, OCCPROB, MINCONSTRAINTDIST, VIS, fig_idx)
if ~exist('VIS', 'var'); VIS = inf; end;
[rows, cols, ~] = size(uv);

% find occd, occr
ind_occ = occ_prob > OCCPROB;
[occd, occr, bad] = find_occlusion_occluder_pairs( ...
  ind_occ, uv, uv_rev, MINCONSTRAINTDIST);
constraint_weights = occ_prob(ind_occ);
constraint_weights(bad) = [];

% find indices:
occd_i = sub2ind([rows, cols], occd(:,1), occd(:,2) );
occr_i = sub2ind([rows, cols], occr(:,1), occr(:,2) );

constraints = [occr_i, occd_i];

if VIS < 150 && fig_idx
  uv_img = im2double(flowToColor(uv));
  occ3_msk = repmat(ind_occ, [1, 1, 3]);
  occ_drawn = sc(double(ind_occ) .* occ_prob, 'jet');
  occ_img = 0.1*uv_img + 0.9*(uv_img.*~occ3_msk + occ_drawn.*occ3_msk);
  cue_img = vis_cues(uv_img, constraints, constraint_weights, 0.0);
  fig(fig_idx); clf; imagesc([occ_img; cue_img]);
end
end
