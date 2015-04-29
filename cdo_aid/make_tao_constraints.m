% returns simple constraints / constraint weights
function [constraints, constraint_weights, indbad] = make_tao_constraints( ...
  occ_prob, uv, uv_rev, OCCPROB, MINCONSTRAINTDIST, VIS, lol)
if ~exist('VIS', 'var'); VIS = inf; end;

[rows, cols, ~] = size(uv);

% find occd, occr
ind_occ = occ_prob > OCCPROB;
[occd, occr, indbad] = find_occlusion_occluder_pairs( ...
  ind_occ, uv, uv_rev, MINCONSTRAINTDIST);
constraint_weights = occ_prob(ind_occ);
constraint_weights(indbad) = [];

% find indices:
occd_i = sub2ind([rows, cols], occd(:,1), occd(:,2) );
occr_i = sub2ind([rows, cols], occr(:,1), occr(:,2) );

constraints = [occr_i, occd_i];

if VIS < 150 && exist('lol', 'var');
  uv_img = im2double(flowToColor(uv));
  occ3_msk = repmat(ind_occ, [1, 1, 3]);
  occ_drawn = sc(double(ind_occ) .* occ_prob, 'jet');
  occ_img = 0.1*uv_img + 0.9*(uv_img.*~occ3_msk + occ_drawn.*occ3_msk);
  % occ_img = uv_img * 0.1 + occ_drawn * 0.9;
  cue_img = vis_cues(uv_img, constraints, constraint_weights, 0.0);
  fig(lol); clf; imagesc([occ_img; cue_img]);
end
end