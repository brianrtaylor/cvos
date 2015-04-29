%--------------------------------------------------------------------
% make_tao_constraint_weights
%
% @param: opts : contains useful options
%   * DO_CONS_PERTURB : 1 to perturb constraints with gmm
%   * CONS_PERTURB : parameters associated with this process
%   * CONSTRAINTS_DIVWEIGHT : 1 to use negative divergence for weighting constraints
%   * relative_weights_constraints : weights between intensity and flow and gpb
%--------------------------------------------------------------------
function [constraints, constraint_weights, w] = make_tao_constraint_weights( ...
  constraints_in, constraint_weights_in, uvb, uvf, I, I_weights, E, opts, lol)
v2struct(opts);

% input
[rows, cols, ~] = size(I);
imsize = [rows, cols];

if opts.DO_CONS_NOW_PERTURB && size(constraints_in, 1) > 0;
  groups = utils_group_constraints(constraints_in, imsize, CONS_PERTURB.GROUP_SIZE);
  fprintf('Computing constraint weights...\n');

  fprintf('Learning GMM.\n');
  tic;
  [local_gmm_bg, local_gmm_fg] = learn_constraint_gmm( ...
    constraints_in, groups, I, max(0.0, 1.0 - I_weights), imsize, opts.CONS_PERTURB);
  toc;
  
  fprintf('Perturbing constraints.\n');
  tic;
  [constraints, gmm_weights] = perturb_constraints(...
    constraints_in, groups, local_gmm_fg, local_gmm_bg, ...
    I, max(0.0, 1.0 - I_weights), CONS_PERTURB);  
  toc;
else
  gmm_weights = ones(size(constraints_in,1),1);
  constraints = constraints_in;
end

% check for constraints falling in the same location, if so, remove them
ind_fail = constraints(:,1) == constraints(:,2);
constraints(ind_fail, :) = [];
constraint_weights_in(ind_fail) = [];
gmm_weights(ind_fail) = [];

if isempty(constraints);
  constraints = [];
  constraint_weights = [];
  w = [];
  
  if VIS < 250 && exist('lol', 'var');
    fig(lol); clf; title('no constraints :(');
  end
else
  
  % incorporate color, flow, gpb edge crossing, divergence into weight
  [w, w_div, ~, w_uvb, w_uvf, w_edgecross] = ...
    make_pixelwise_nonadj_weights(constraints, ...
    I, uvb, uvf, E, relative_weights_constraints, true); % used to be Ilab
  if DO_CONS_NOW_PERTURB;
    %  constraint_weights = constraint_weights .* gmm_weights;
    b = relative_weights_constraints / sum(relative_weights_constraints);
    w = b(1) * (1.0 - gmm_weights) ...
      + b(2) * ((w_uvb + w_uvf) / 2.0) ...
      + b(4) * w_edgecross;
  end
  w = 1.0 - w;
  constraint_weights = w .* constraint_weights_in;
  if DO_CONSTRAINT_DIVWEIGHT;
    constraint_weights = constraint_weights .* w_div;
  end
  
  
  if VIS < 250 && exist('lol', 'var');
    uvf_img  = im2double(flowToColor(uvf));
    cue_img1 = vis_cues(uvf_img, constraints_in, constraint_weights_in, 0.0);
    cue_img2 = vis_cues(uvf_img, constraints, constraint_weights, 0.0);
    fig(lol); clf; imagesc([cue_img1; cue_img2]);
  end
end
end
