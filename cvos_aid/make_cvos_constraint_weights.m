%-----------------------------------------------------------------------------
% make_cvos_constraints weights
%
% computes the weights for the occluder-occluded constraints based on motion,
% optical flow, and edge strength 
%
% @return: cons: output occluder-occluded constraints
% @return: cons_weights: output weights
% @return: w: modifier for existing weights
% @param: cons_in: list of occluder-occluded constraints
% @param: cons_weights_in: initial weights computed from occ strength
% @param: uvb: backward warping from current frame to the frame before
% @param: uvf: forward warping from current frame to the frame after
% @param: I: current image frame
% @param: I_weights: weights computed on current image frame
% @param: E: edge strength
% @param: opts: struct of useful options
%   * DO_CONS_PERTURB: 1 to perturb constraints with gmm
%   * CONS_PERTURB: parameters associated with this process
%   * CONSTRAINTS_DIVWEIGHT: to use negative divergence to weight constraints
%   * relative_weights_constraints: weights between intensity and flow and gpb
% @param: lol: figure number for visualization
%-----------------------------------------------------------------------------
function [cons, cons_weights, w] = make_cvos_constraint_weights(cons_in, ...
  cons_weights_in, uvb, uvf, I, I_weights, E, opts, lol)
v2struct(opts);
[rows, cols, ~] = size(I); imsize = [rows, cols];
if opts.DO_CONS_NOW_PERTURB && size(cons_in, 1) > 0;
  fprintf('Computing constraint weights\n');
  groups = utils_group_constraints(cons_in, imsize, CONS_PERTURB.GROUP_SIZE);
  fprintf('Learning GMM\n');
  [local_gmm_bg, local_gmm_fg] = learn_constraint_gmm( ...
    cons_in, groups, I, max(0.0, 1.0 - I_weights), imsize, opts.CONS_PERTURB);
  fprintf('Perturbing constraints\n');
  [cons, gmm_weights] = perturb_constraints(cons_in, groups, local_gmm_fg, ...
    local_gmm_bg, I, max(0.0, 1.0 - I_weights), CONS_PERTURB);  
else
  gmm_weights = ones(size(cons_in,1),1);
  cons = cons_in;
end

% check for constraints falling in the same location, if so, remove them
ind_fail = cons(:,1) == cons(:,2);
cons(ind_fail, :) = [];
cons_weights_in(ind_fail) = [];
gmm_weights(ind_fail) = [];

if isempty(cons);
  cons = [];
  cons_weights = [];
  w = [];
  
  if VIS < 250 && exist('lol', 'var');
    fig(lol); clf; title('no constraints :(');
  end
else
  % incorporate color, flow, gpb edge crossing, divergence into weight
  [w, w_div, ~, w_uvb, w_uvf, w_edgecross] = ...
    make_pixelwise_nonadj_weights(cons, I, uvb, uvf, E, ...
    relative_weights_constraints, true);
  if DO_CONS_NOW_PERTURB;
    b = relative_weights_constraints / sum(relative_weights_constraints);
    w = b(1) * (1.0 - gmm_weights) ...
      + b(2) * ((w_uvb + w_uvf) / 2.0) ...
      + b(4) * w_edgecross;
  end
  w = 1.0 - w;
  cons_weights = w .* cons_weights_in;
  if DO_CONSTRAINT_DIVWEIGHT;
    cons_weights = cons_weights .* w_div;
  end
  
  if VIS < 250 && exist('lol', 'var');
    uvf_img  = im2double(flowToColor(uvf));
    cue_img1 = vis_cues(uvf_img, cons_in, cons_weights_in, 0.0);
    cue_img2 = vis_cues(uvf_img, cons, cons_weights, 0.0);
    fig(lol); clf; imagesc([cue_img1; cue_img2]);
  end
end
end
