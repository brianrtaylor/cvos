% returns simple constraints / constraint weights
function [constraints_uno, constraint_weights_uno, indbad_uno, constraints_uno_img] = make_cvos_unary_constraints( ...
  occ_prob, uv, OCCPROB)

[rows, cols, ~] = size(uv);

% find pixels now occd, future occr
ind_occ = occ_prob > OCCPROB;
[occr_uno, indbad_uno, occr_uno_img] = find_occlusion_future_occluder_uno( ...
  ind_occ, uv);
constraint_weights_uno = occ_prob(ind_occ);
constraint_weights_uno(indbad_uno) = [];

% find indices:
occr_uno_i = sub2ind([rows, cols], occr_uno(:,1), occr_uno(:,2));

constraints_uno = occr_uno_i;

[constraints_uno, constraint_weights_uno, ind_uno] = aggregate_list( ...
  constraints_uno, constraint_weights_uno, false);

constraints_uno_img = zeros(rows, cols);
constraints_uno_img(constraints_uno) = constraint_weights_uno;
end
