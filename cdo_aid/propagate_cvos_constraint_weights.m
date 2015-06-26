%--------------------------------------------------------------------
% propagate_cvos_constraint_weights
%
% @param: opts : options and parameters
%   * DO_CONS_PERTURB : 1 to perturb constraints
%   * CONSTRAINT_INSIDE_FG_MEMORY_FACTOR : falloff factor for constraints 
%       inside of a segmentation
%   * DO_CONSTRAINT_WARP_WEIGHT : TODO: change to DO_C..., 1 if we penalize 
%       constraint weight for warping
%   * CONSTRAINT_MEMORY_FACTOR : falloff factor of constraints everywhere
%   * CONSTRAINT_WEIGHT_THRESH : minimum weight for constraints. those 
%       below this are elminated
%   * DO_RECOMPUTE_OLD_WEIGHTS_IN_CURRENT_FRAME : recompute pixel-wise 
%       weights in the current frame
%   * CONSTRAINT_DIVWEIGHT : 
%   * relative_constraint_weights : proportion of weight from intensity, 
%       motion difference, gpb, etc.
%   * GMMPKG : parameters and settings for running gmm
%   * CAUSAL : 1 if causal setup (not frame by frame)
%--------------------------------------------------------------------
function [constraints, constraint_weights, constraints_causal, constraint_weights_old, valid] = ...
  propagate_cvos_constraint_weights( ...
  I0, I1, i1, uvb, uvf, past, past_constraints, ...
  past_constraint_weights, prev_weights, weights, ...
  xi, E, constraints_now, constraint_weights_now, opts, objects, lol)

v2struct(opts);
[rows, cols, ~] = size(I1);
imsize = [rows, cols];

%------------------------------------------------------------------
% meat
%------------------------------------------------------------------
constraints_causal = [];
constraint_weights_old = past_constraint_weights;
valid = ones(size(past_constraints, 1), 1);

%------------------------------------------------------------------
% stop early if not causal or no constraints
%------------------------------------------------------------------
if isempty(past_constraints) || ~CAUSAL;
  constraints = constraints_now;
  constraint_weights = constraint_weights_now;
  return;
end
  
%------------------------------------------------------------------
% gmm warping of old constraints
%------------------------------------------------------------------   
if DO_CONS_PERTURB;
  [constraints_causal, gmm_weights, warp_weights, valid] = ...
    gmm_constraint_warping(past.uvf_cbf_bflt, I0, I1, ...
    prev_weights, weights, past_constraints, ...
    past.layers, GMMPKG, WARP_SAFESPEEDSQUARED);
else
  [constraints_causal, warp_weights, valid] = ...
    warp_constraints_forward(past_constraints, past.uvf_cbf_bflt, ...
    WARP_SAFESPEEDSQUARED); % uvf_cbf or uvf
end
constraint_weights_old(~valid) = [];

%------------------------------------------------------------------
% stop early if no valid constraints
%------------------------------------------------------------------
if sum(valid) <= 0;
  constraints = [];
  constraint_weights = [];
  return;
end
  
%------------------------------------------------------------------
% diminish constraints by age and warping
%------------------------------------------------------------------
w_constraints = 1;

% 1. penalty for moving (larger motion = greater penalty)
if DO_CONSTRAINT_WARP_WEIGHT;
  w_constraints = w_constraints .* warp_weights;
end

% 2. computing penalty for motion different from object motion
nobjects = size(objects, 1);
if DO_CONSTRAINT_OBJECT_WARP_WEIGHT && nobjects > 0;
  ind_occr = past_constraints(:, 1);
  w_fg_occr_uvf = past.w_fg_uvf(ind_occr);
    
  if VIS < 120;
    w_constraints2 = w_constraints .* w_fg_occr_uvf(valid);
    vcn2 = vis_cues(i1, constraints_causal, w_constraints2);
    vcn4 = vis_cues(i1, constraints_causal, w_constraints);
    fig(lol); clf; imagesc([vcn2;vcn4]);
  end

  w_constraints = w_constraints .* w_fg_occr_uvf(valid);
end

% age discount
w_constraints = (1 - CONSTRAINT_MEMORY_FACTOR) .* w_constraints;

%------------------------------------------------------------------
% boost constraints that agree with segmentation
%------------------------------------------------------------------  
% compute the weight of the good constraints
[w, ~, ~, w_uvb, w_uvf, w_edgecross] = ...
  make_pixelwise_nonadj_weights(constraints_causal, ...
  I1, uvb, uvf, E, relative_weights_constraints, false);
if DO_CONS_PERTURB;
  b = relative_weights_constraints / sum(relative_weights_constraints);
  w = b(1) * (1.0 - gmm_weights) ...
    + b(2) * ((w_uvb + w_uvf) / 2.0) ...
    + b(4) * w_edgecross;
end
good_weights = 1.0 - w;

if DO_RECOMPUTE_OLD_WEIGHTS_IN_CURRENT_FRAME;
  constraint_weights_old = good_weights;
end

% add segmentation bonus
xi(~valid) = [];
good_constraints = max(0.0, 1 - xi) .* good_weights;
bad_constraints = (max(0.0, xi) .* good_weights) / 2;

% finally actually multiplies the constraint weights out
constraint_weights_old2 = max(0.0, constraint_weights_old .* w_constraints ...
  + good_constraints - bad_constraints);
  
if VIS < 220; 
  vcn1 = vis_cues(i1, constraints_causal, min(1.0, constraint_weights_old));
  vcn2 = vis_cues(i1, constraints_causal, w_constraints);
  vcn3 = vis_cues(i1, constraints_causal, good_constraints);
  vcn4 = vis_cues(i1, constraints_causal, min(1.0, constraint_weights_old2));
  fig(100 + lol); clf; imagesc([[vcn1, vcn2];[ vcn3, vcn4]]);
end

constraint_weights_old = constraint_weights_old2;

%------------------------------------------------------------------
% final aggregation
%------------------------------------------------------------------
[constraints, constraint_weights, ~] = aggregate_pairs_fast( ...
  double([constraints_now; constraints_causal]), ...
  double([constraint_weights_now; constraint_weights_old]));
  
%------------------------------------------------------------------
% visualizations
%------------------------------------------------------------------
if VIS < 180 && sum(valid) > 0;
  a00 = vis_weights_occlusions(past.weights, past_constraints, imsize, ...
    past_constraint_weights);
  aaa = vis_weights_occlusions(weights, constraints_causal, imsize, ...
    constraint_weights_old);
  bbb = vis_weights_occlusions(weights, constraints_now, imsize, ...
    constraint_weights_now);
  ccc = vis_weights_occlusions(weights, constraints_causal, imsize, ...
    warp_weights);
  ddd = vis_weights_occlusions(weights, constraints, imsize, ...
    constraint_weights);
  fig(410); clf; imagesc([a00, aaa, bbb; aaa * 0, ccc, ddd]); axt;
end

%------------------------------------------------------------------
% some weight based pruning on constraints
%------------------------------------------------------------------
weak_inds = constraint_weights < CONSTRAINT_WEIGHT_THRESH;
constraints(weak_inds, :) = [];
constraint_weights(weak_inds) = [];
end
