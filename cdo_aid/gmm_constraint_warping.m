% -------------------------------------------------------------------------
% INPUT:
% -------------------------------------------------------------------------
%
%   uvf_prev = uvf(t-1)
%   Ib = I(t-1) (previous image)
%   I0 = I(t) (current image)
%   prev_weights = weights(t-1) (weights in previous frame)
%   weights = weights(t) (current weights)
%   constraints_prev - constraints defined in the domain of the previous
%       frame.
%   prev_layers = segmentation(t-1)
%   params - parameter structure.............
%
% -------------------------------------------------------------------------
% OUTPUT:
% -------------------------------------------------------------------------
%
%   constraints - in the domain of the current frame, appropriately
%   deformed...
%   gmm_weights - weights, computed according to:
%           w(i) = 0.5*( Pr(x_{FG} \in FG) + Pr(x_{BG} \in BG ) )
%       if weight is near 1, the constraint is reliable.
%       if weight is near 0.5, the constraint is pretty unreliable.
%       it is fairly unlikely that weight reaches low values (~0.1 and below),
%       but if it does, background GMM is really bad at modeling background
%       pixel intensity, and foreground GMM is really bad at modeling
%       foreground pixel...
%   warp_weights - weights
%   valid - which of the input constraints / weights are valid (1) or not
%       (0), in order as input
%   ind_out - output ordering of input constraints. so constraints =
%       constraints_prev(ind_out, :), essentially, assuming none invalid
%
function [constraints, gmm_weights, warp_weights, valid] = ...
    gmm_constraint_warping( ...
  uvf_prev, Ib, I0, prev_weights, weights, constraints_prev, ...
  prev_layers, params, WARP_SAFESPEEDSQUARED)

    N = size(constraints_prev, 1);
    valid = false(N, 1); % fails through, all fail
    % ---------------------------------------------------------------------
    %       Block that adds previous constraints to the current ones
    % ---------------------------------------------------------------------
    imsize = [size(I0,1), size(I0,2)];

    constraint_prev_vals = prev_layers( constraints_prev(:,1) ) ...
      - prev_layers( constraints_prev(:,2) );
    
    % >= 1 --> constraint was satisfied.
    % == 0 --> unsatisfied (this cost something)
    %  < 0 --> constraint was determined to be incorrect

    % -----------------------------------------------------------------
    %             Split constraints into two groups
    %   1. active constraints: c_{t-1}(y) > c_{t-1}(x)  
    %   1. inactive constraints: c_{t-1}(y) = c_{t-1}(x)        
    % -----------------------------------------------------------------
    tokeep = (constraint_prev_vals >= 0.5);
    valid_active = zeros(sum(tokeep), 1);
    ind_active = find(tokeep);
    if sum(tokeep) > 0
        constraints_prev_active = constraints_prev(tokeep,:);
        if params.cons_perturb.DO_CONS_PERTURB;

          groups = utils_group_constraints(constraints_prev_active, imsize, params.cons_perturb.GROUP_SIZE);

          fprintf('Learning GMM.\n');
          tic;
          [local_gmm_bg_active, local_gmm_fg_active] = ...
            learn_constraint_gmm(constraints_prev_active, groups, ...
            Ib, max(0.0, 1.0 - prev_weights), imsize, params.cons_perturb, ...
            prev_layers);
          toc;
          [constraints_prev_active, warp_weights_active, valid_active] = ...
            warp_constraints_forward(constraints_prev_active, uvf_prev, ...
            WARP_SAFESPEEDSQUARED);
          % Validate after warping:
          groups = groups(valid_active);

          fprintf('Perturbing constraints.\n');
          tic;
          [constraints_prev_active, gmm_weights_active] = ...
            perturb_constraints(constraints_prev_active, groups, ...
            local_gmm_fg_active, local_gmm_bg_active, I0, ...
            max(0.0, 1.0 - weights), params.cons_perturb);
          toc;          
        else
            [constraints_prev_active, warp_weights_active, valid_active] = ...
              warp_constraints_forward(constraints_prev_active, uvf_prev, ...
              WARP_SAFESPEEDSQUARED);                     
            gmm_weights_active = ones( size(constraints_prev_active, 1), 1);
        end
        valid(ind_active) = valid_active;
    else
      constraints_prev_active = zeros(0,2);
      warp_weights_active = [];
      gmm_weights_active = [];
    end

    tokeep = (constraint_prev_vals >= -0.5 & constraint_prev_vals < 0.5);
    valid_inactive = zeros(sum(tokeep), 1);
    ind_inactive = find(tokeep);
    if sum(tokeep) > 0
        constraints_prev_inactive = constraints_prev(tokeep,:);
        if params.cons_perturb.DO_CONS_PERTURB;

          groups = utils_group_constraints(constraints_prev_inactive, imsize, params.cons_perturb.GROUP_SIZE);

          fprintf('Learning GMM.\n');
          tic;
          [local_gmm_bg_inactive, local_gmm_fg_inactive] = ...
            learn_constraint_gmm(constraints_prev_inactive, groups, ...
            Ib, max(0.0, 1.0 - prev_weights), imsize, params.cons_perturb);
          toc;
          [constraints_prev_inactive, warp_weights_inactive, valid_inactive] = ...
            warp_constraints_forward(constraints_prev_inactive, uvf_prev, ...
            WARP_SAFESPEEDSQUARED);
          groups = groups(valid_inactive);

          fprintf('Perturbing constraints.\n');
          tic;
          [constraints_prev_inactive, gmm_weights_inactive] = ...
            perturb_constraints(constraints_prev_inactive, groups, ...
            local_gmm_fg_inactive, local_gmm_bg_inactive, I0, ...
            max(0.0, 1.0 - weights), params.cons_perturb);
          toc;
        else
          [constraints_prev_inactive, warp_weights_inactive, valid_inactive] = ...
            warp_constraints_forward(constraints_prev_inactive, uvf_prev, ...
            WARP_SAFESPEEDSQUARED);
          if isempty(valid); valid = valid_inactive; else valid = valid & valid_inactive; end;
          gmm_weights_inactive = ones( size(constraints_prev_inactive, 1), 1);                
        end
        valid(ind_inactive) = valid_inactive;
    else
        constraints_prev_inactive = zeros(0,2);
        warp_weights_inactive = [];
        gmm_weights_inactive = [];
    end
    fprintf('# constraints: active: %d inactive: %d\n', ...
        size(constraints_prev_active,1), size(constraints_prev_inactive,1) );
    %------------------------------------------------------------------
    %                Combine two groups of constraints
    %------------------------------------------------------------------
    constraints = zeros(N, 2);
    constraints(ind_inactive(valid_inactive), :) = constraints_prev_inactive;
    constraints(ind_active(valid_active), :) = constraints_prev_active;
    constraints = constraints(logical(valid), :);
    
    gmm_weights = zeros(N, 1);
    gmm_weights(ind_inactive(valid_inactive)) = gmm_weights_inactive;
    gmm_weights(ind_active(valid_active)) = gmm_weights_active;
    gmm_weights = gmm_weights(logical(valid), :);
    gmm_weights(isnan(gmm_weights) | isinf(gmm_weights)) = 0.0;
    
    warp_weights = zeros(N, 1);
    warp_weights(ind_inactive(valid_inactive)) = warp_weights_inactive;
    warp_weights(ind_active(valid_active)) = warp_weights_active;
    warp_weights = warp_weights(logical(valid), :);
    warp_weights(isnan(warp_weights) | isinf(warp_weights)) = 0.0;
end
