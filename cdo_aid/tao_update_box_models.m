function boxes = tao_update_box_models( ...
  boxes, layers, i1_bflt, Dx, Dy, weights, opts)
BOX_LIMIT = inf;
GMM_PROB_FG_MASK_THRESH = 0.75;
GMM_PROB_BG_MASK_THRESH = 0.25;
FGTHRESH = 0.75;

CP = opts.CONS_PERTURB;

if opts.CAUSAL && opts.BOXHELP;
  tic
  
  % TODO: jitter box towards layer ege boundary
  [new_boxes, boxes, valid_boxes] = bboxes_for_mask_update( ...
    layers, Dx, Dy, opts.BOX_RAD, boxes, BOX_LIMIT);
  n_newboxes = length(new_boxes);
  n_boxes = length(boxes);
  all_boxes = [boxes; new_boxes];
  
  toc
  
  %-----------------------------------------------------------------
  % learn gmm if foreground objects exist and we have valid boxes
  %-----------------------------------------------------------------
  if ~isempty(all_boxes) && max(layers(:)) >= 1;

    [mu_fg, cov_fg, pi_fg, mu_bg, cov_bg, pi_bg] = learn_bbox_gmm_mex( ...
      i1_bflt, double(layers), double(weights), all_boxes, ...
      CP.NUM_GMM_CLUSTERS, CP.NUM_GMM_REPETITIONS, CP.NUM_GMM_ITERATIONS, ...
      GMM_PROB_FG_MASK_THRESH, GMM_PROB_BG_MASK_THRESH);
    
    %-----------------------------------------------------------------
    % choose which gmm to use (model selection on fg size)
    %-----------------------------------------------------------------
    % for new boxes, learn fresh
    if n_newboxes > 0;
      % if isempty(all_boxes(1).fg_gmm_mu); % hasn't been set yet
      for bk = 1:n_newboxes;
        new_boxes(bk).gmm_change = true;
        new_boxes(bk).fg_gmm_mu = mu_fg(:,:,bk);
        new_boxes(bk).fg_gmm_cov = cov_fg(:,:,bk);
        new_boxes(bk).fg_gmm_pi = pi_fg(:,bk);
        new_boxes(bk).bg_gmm_mu = mu_bg(:,:,bk);
        new_boxes(bk).bg_gmm_cov = cov_bg(:,:,bk);
        new_boxes(bk).bg_gmm_pi = pi_bg(:,bk);
      end
      
      [box_fg, box_bg] = eval_gmm_bboxes_mex(i1_bflt, new_boxes);
      for bk = 1:n_newboxes;
        new_boxes(bk).fg_prob = box_fg(bk).prob;
        new_boxes(bk).bg_prob = box_bg(bk).prob;
      end
    end
      
    if n_boxes > 0; % for old boxes, do some checks
      % evaluate historic model (frame 0 model)
      for bk = 1:n_boxes;
        boxes(bk).fg_gmm_mu(:, (CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).fg_gmm_cov(:, (CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).fg_gmm_pi((CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        
        boxes(bk).bg_gmm_mu(:, (CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).bg_gmm_cov(:, (CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).bg_gmm_pi((CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
      end
      [box_fg_1, box_bg_1] = eval_gmm_bboxes_mex(i1_bflt, boxes);
        
      % create and evaluate updated model
      up_boxes = boxes;
      for bk = 1:n_boxes;
        up_boxes(bk).fg_gmm_mu  = [boxes(bk).fg_gmm_mu, mu_fg(:,:,bk)];
        up_boxes(bk).fg_gmm_cov = [boxes(bk).fg_gmm_cov, cov_fg(:,:,bk)];
        up_boxes(bk).fg_gmm_pi  = [boxes(bk).fg_gmm_pi; pi_fg(:,bk)];
        
        up_boxes(bk).bg_gmm_mu  = [boxes(bk).bg_gmm_mu, mu_bg(:,:,bk)];
        up_boxes(bk).bg_gmm_cov = [boxes(bk).bg_gmm_cov, cov_bg(:,:,bk)];
        up_boxes(bk).bg_gmm_pi  = [boxes(bk).bg_gmm_pi; pi_bg(:,bk)];
      end
      [box_fg_2, box_bg_2] = eval_gmm_bboxes_mex(i1_bflt, up_boxes);
  
      % make selection of historic vs updated model
      for bk = 1:n_boxes;
        fg1 = box_fg_1(bk).prob > FGTHRESH;
        fg2 = box_fg_2(bk).prob > FGTHRESH;
        
        nfg1 = sum(fg1(:));
        nfg2 = sum(fg2(:));
        
        % TODO: maybe keep whichever segmentation is smoother +
        % smaller, define some energy to evaluate it
        % classic model has larger fg, then keep updated region
        if (nfg1 > nfg2);
          boxes(bk).gmm_change = true;
          
          boxes(bk).fg_gmm_mu  = up_boxes(bk).fg_gmm_mu;
          boxes(bk).fg_gmm_cov = up_boxes(bk).fg_gmm_cov;
          boxes(bk).fg_gmm_pi  = up_boxes(bk).fg_gmm_pi;
          
          boxes(bk).bg_gmm_mu  = up_boxes(bk).bg_gmm_mu;
          boxes(bk).bg_gmm_cov = up_boxes(bk).bg_gmm_cov;
          boxes(bk).bg_gmm_pi  = up_boxes(bk).bg_gmm_pi;
          
          boxes(bk).fg_prob = box_fg_2(bk).prob;
          boxes(bk).bg_prob = box_bg_2(bk).prob;
          
          % update patches for colour model update
          % I think the conf requires both the layers patch at the time
          % and the confidence, so it can be recomputed (TODO)
          % all_boxes(bk).patches = cat(3, all_boxes(bk).patches, ...
          %   i1_bflt(ys, xs, :));
        else
          boxes(bk).gmm_change = false;
          boxes(bk).fg_prob = box_fg_1(bk).prob;
          boxes(bk).bg_prob = box_bg_1(bk).prob;
        end
      end
    end
      
    all_boxes = [boxes; new_boxes];
    
    %-------------------------------
    % prune out any bad NaN gmms
    %-------------------------------
    % after choosing which gmms
    bad_gmm_bg = isnan(sum(pi_bg, 1));
    bad_gmm_fg = isnan(sum(pi_fg, 1));
    
    % assert(bad_gmm_bg == bad_gmm_fg, 'the same box is bad');
    bad_gmm = bad_gmm_fg | bad_gmm_bg;
    
    all_boxes(bad_gmm) = [];
    
    % % % %-------------------------------
    % % % % learn the flow model
    % % % %-------------------------------
    % % % % learning flow model (mean fg flow, mean bg flow)
    % % % [uv_bg, conf_uv_bg, uv_fg, conf_uv_fg, valid_uv] = learn_bbox_uv_model( ...
    % % %   all_boxes, uvf_bflt);
    % % % all_boxes = all_boxes(valid_uv);
    % % % n_all_boxes = size(all_boxes, 1);
    % % % for bk = 1:n_all_boxes;
    % % %   all_boxes(bk).u_bg = uv_bg(bk, 1);
    % % %   all_boxes(bk).v_bg = uv_bg(bk, 2);
    % % %   all_boxes(bk).u_fg = uv_fg(bk, 1);
    % % %   all_boxes(bk).v_fg = uv_fg(bk, 2);
    % % %   all_boxes(bk).conf_uv_bg = conf_uv_bg(bk); % some variance
    % % %   all_boxes(bk).conf_uv_fg = conf_uv_fg(bk); % some variance
    % % % end
    
    %-----------------------------------------------------------------
    % evaluate on current frame to obtain confidences for new boxes
    %-----------------------------------------------------------------
    [conf_colour_now, conf_colour_denom] ...
      = compute_colour_model_bbox_conf(all_boxes, layers);
      
    for bk = 1:length(all_boxes);
      if isempty(all_boxes(bk).conf_colour);
        all_boxes(bk).conf_colour_h = conf_colour_now(bk);
        all_boxes(bk).conf_colour_denom_h = conf_colour_denom(bk);
      end
      all_boxes(bk).conf_colour = conf_colour_now(bk);
      all_boxes(bk).conf_colour_denom = conf_colour_denom(bk);
    end
    
    % finish it
    boxes = all_boxes;
  else
    % TODO: do somerthing smarter, or make this a parameter
    % box_conf = box_conf * 0.9;
  end
end
end
