%--------------------------------------------------------------------
% propagate_tao_constraint_weights
%
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
function [constraints, constraint_weights, constraints_causal, constraint_weights_old, valid] = propagate_tao_constraint_weights( ...
  I0, I1, i1, uvb, uvf, past, past_constraints, ...
  past_constraint_weights, prev_weights, weights, ...
  xi, E, constraints_now, constraint_weights_now, opts, ...
  prob_box_fg, count_box_fg, objects, lol)

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
% vasiliy's gmm warping of old constraints
%------------------------------------------------------------------   
if DO_CONS_PERTURB;
  fake = ones(size(past_constraints, 1), 1);
  [constraints_causal, ~, gmm_weights, warp_weights, valid] = ...
    gmm_constraint_warping(past.uvf_cbf_bflt, I0, I1, ...
    prev_weights, weights, past_constraints, ...
    fake, past.layers, GMMPKG, WARP_SAFESPEEDSQUARED);
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

% % % gmm confidence discount % (using this in weights instead)
% % if DO_CONS_PERTURB;
% %   w_constraints = w_constraints .* gmm_weights;
% % end

% age discount
w_constraints = (1 - CONSTRAINT_MEMORY_FACTOR) .* w_constraints;

% % % %------------------------------------------------------------------
% % % % constraints refinement (temporally consistent layer validation)
% % % %------------------------------------------------------------------
% % %   % TODO: ISN'T CORRECT, when layers switch from one frame to the next,
% % %   % like when ppl overlap (people2)
% % %   % TODO: make it work more for when an object occludes another, less so
% % %   % just the "falling off" thing
% % %   % TODO: can we assume. this is okay, because when an object switches
% % %   % layers completely, then there should be enough occlusions from that
% % %   % frame itself to create all the necessary things for this object to pop.
% % %   % we don't need to have historical evidence
% % %   LAYER_T0_VALIDATION = true;
% % %   % if the layers don't stay the same after warping, assume constraint
% % %   % "fell-off" the object and doesn't belong
% % %   lay_occr = past.layers(past_constraints(:, 1));
% % %   lay_occr(~valid) = [];
% % %   if LAYER_T0_VALIDATION;
% % %     layer_t0_mask = past.t0.layers;
% % %     lay_t0_occr = layer_t0_mask(constraints_causal(:, 1));
% % %     lay_t0_occd = layer_t0_mask(constraints_causal(:, 2));
% % %     w_valid = max(0.0, 1.0 - abs(lay_t0_occr - lay_occr));
% % %
% % %     % layer_t0_mask_g = imfilter(layer_t0_mask, g_tiny);
% % %     layer_t0_mask_g = imdilate(layer_t0_mask, dsk);
% % %     lay_t0_occr_g = layer_t0_mask_g(constraints_causal(:, 1));
% % %     lay_t0_occd_g = layer_t0_mask_g(constraints_causal(:, 2));
% % %     w_valid_g = max(0.0, 1.0 - abs(lay_t0_occr_g - lay_occr));
% % %
% % %     % TODO: look into this. the visualization shows it seems to
% % %     % not be working, but I think this should kill off some of
% % %     % those occlusions trailing behind the car rear top window
% % %     % in cars 1 frame 3 or so. w_valid should be true for those
% % %     % constraints that stay on the object, and should be 0 for
% % %     % those constraints that "fall off" the object from the
% % %     % warping. And those behind the window should "fall off"
% % %     w_valid = max(w_valid, w_valid_g);
% % %
% % %     w_constraints = w_constraints .* w_valid;
% % %   end
  
%------------------------------------------------------------------
% boost constraints that agree with segmentation
%------------------------------------------------------------------  
% compute the weight of the good constraints
[w, ~, ~, w_uvb, w_uvf, w_edgecross] = ...
  make_pixelwise_nonadj_weights(constraints_causal, ...
  I1, uvb, uvf, E, relative_weights_constraints, false);
if DO_CONS_PERTURB;
  %  constraint_weights = constraint_weights .* gmm_weights;
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
  
% some quick visuals
if VIS < 220; 
  vcn1 = vis_cues(i1, constraints_causal, min(1.0, constraint_weights_old));
  vcn2 = vis_cues(i1, constraints_causal, w_constraints);
  vcn3 = vis_cues(i1, constraints_causal, good_constraints);
  vcn4 = vis_cues(i1, constraints_causal, min(1.0, constraint_weights_old2));
  fig(100 + lol); clf; imagesc([[vcn1, vcn2];[ vcn3, vcn4]]);
end

constraint_weights_old = constraint_weights_old2;

% %------------------------------------------------------------------
% % attempt to remove constraints from self-disocclusion region
% %------------------------------------------------------------------ 
%   % if both have high foreground prob in current frame, likely it's a
%   % self-occlusion / disocclusion, and both belong to foreground
%   if exist('prob_box_fg', 'var') && exist('count_box_fg', 'var');
%     count_box_fg_occr = count_box_fg(constraints_causal(:, 1));
%     count_box_fg_occd = count_box_fg(constraints_causal(:, 2));
%
%     prob_box_fg_occr = prob_box_fg(constraints_causal(:, 1));
%     prob_box_fg_occd = prob_box_fg(constraints_causal(:, 2));
%
%     % TODO: check this, case we should exp it or something instead
%     prob_box_fg_constraint = max(0.0, prob_box_fg_occr - prob_box_fg_occd);
%     prob_box_fg_constraint((count_box_fg_occr == 0) | (count_box_fg_occd == 0)) = 1.0;
%
%     constraint_weights_old = constraint_weights_old .* prob_box_fg_constraint;
%   end

%------------------------------------------------------------------
% final aggregation
%------------------------------------------------------------------
[constraints, constraint_weights, constraint_inds] = aggregate_pairs_fast( ...
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
%   
%   iii = vis_weights_occlusions(weights, constraints_causal, ...
%     imsize, w_valid);
%   hhh = vis_weights_occlusions(weights, constraints_causal, ...
%     imsize, max(0.0, 1 - xi));
%   
%   fff = vis_weights_occlusions(weights, constraints_causal, ...
%     imsize, lay_t0_occr);
%   ggg = vis_weights_occlusions(past.weights, past_constraints(valid, :), ...
%     imsize, lay_occr(valid));
%   
%   fff = vis_weights_occlusions(weights, constraints_causal, ...
%     imsize, w_valid);
%   ggg = vis_weights_occlusions(past.weights, past_constraints(valid, :), ...
%     imsize, w_valid_g);
%   
%   fff = vis_weights_occlusions(weights, constraints_causal, imsize);
%   ggg = vis_weights_occlusions(past.weights, past_constraints(valid, :), ...
%     imsize);
  fig(410); clf; imagesc([a00, aaa, bbb; aaa * 0, ccc, ddd]); axt;
end

%------------------------------------------------------------------
% some weight based pruning on constraints
%------------------------------------------------------------------
weak_inds = constraint_weights < CONSTRAINT_WEIGHT_THRESH;
constraints(weak_inds, :) = [];
constraint_weights(weak_inds) = [];
end
