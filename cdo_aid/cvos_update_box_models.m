function boxes = cvos_update_box_models( ...
  boxes, layers, i1_bflt, Dx, Dy, weights, opts)
BOX_LIMIT = inf;
GMM_PROB_FG_MASK_THRESH = 0.75;
GMM_PROB_BG_MASK_THRESH = 0.25;
FGTHRESH = 0.75;

CP = opts.CONS_PERTURB;

if opts.CAUSAL && opts.BOXHELP;
  [new_boxes, boxes] = bboxes_for_mask_update( ...
    layers, Dx, Dy, opts.BOX_RAD, boxes, BOX_LIMIT);
  n_newboxes = length(new_boxes);
  n_boxes = length(boxes);
  
  %-----------------------------------------------------------------
  % learn gmm if foreground objects exist and we have valid boxes
  %-----------------------------------------------------------------
  if max(layers(:)) >= 1;       
    %-----------------------------------------------------------------
    % choose which gmm to use (model selection on fg size)
    %-----------------------------------------------------------------
    % for new boxes, learn fresh
    if n_newboxes > 0;
      [mu_fg, cov_fg, pi_fg, mu_bg, cov_bg, pi_bg] = learn_bbox_gmm_mex( ...
        i1_bflt, double(layers), double(weights), new_boxes, ...
        CP.NUM_GMM_CLUSTERS, CP.NUM_GMM_REPETITIONS, CP.NUM_GMM_ITERATIONS, ...
        GMM_PROB_FG_MASK_THRESH, GMM_PROB_BG_MASK_THRESH);
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
      [mu_fg, cov_fg, pi_fg, mu_bg, cov_bg, pi_bg] = learn_bbox_gmm_mex( ...
        i1_bflt, double(layers), double(weights), boxes, ...
        CP.NUM_GMM_CLUSTERS, CP.NUM_GMM_REPETITIONS, CP.NUM_GMM_ITERATIONS, ...
        GMM_PROB_FG_MASK_THRESH, GMM_PROB_BG_MASK_THRESH);
        
      % evaluate historic model (frame 0 model)
      for bk = 1:n_boxes;
        boxes(bk).fg_gmm_mu(:, (CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).fg_gmm_cov(:, (CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).fg_gmm_pi((CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).fg_gmm_pi = boxes(bk).fg_gmm_pi/sum(boxes(bk).fg_gmm_pi);
        
        boxes(bk).bg_gmm_mu(:, (CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).bg_gmm_cov(:, (CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).bg_gmm_pi((CP.NUM_GMM_CLUSTERS + 1):end, :) = [];
        boxes(bk).bg_gmm_pi = boxes(bk).bg_gmm_pi/sum(boxes(bk).bg_gmm_pi);
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

        valid = ~isnan(up_boxes(bk).fg_gmm_pi) & ... 
                ~isnan(up_boxes(bk).bg_gmm_pi);
        up_boxes(bk).fg_gmm_mu  = up_boxes(bk).fg_gmm_mu(:, valid);
        up_boxes(bk).fg_gmm_cov = up_boxes(bk).fg_gmm_cov(:, valid);
        up_boxes(bk).fg_gmm_pi  = up_boxes(bk).fg_gmm_pi(valid);        
        up_boxes(bk).fg_gmm_pi = up_boxes(bk).fg_gmm_pi/sum(up_boxes(bk).fg_gmm_pi);

        up_boxes(bk).bg_gmm_mu  = up_boxes(bk).bg_gmm_mu(:, valid);
        up_boxes(bk).bg_gmm_cov = up_boxes(bk).bg_gmm_cov(:, valid);
        up_boxes(bk).bg_gmm_pi  = up_boxes(bk).bg_gmm_pi(valid);        
        up_boxes(bk).bg_gmm_pi = up_boxes(bk).bg_gmm_pi/sum(up_boxes(bk).bg_gmm_pi);
        
        if isempty(up_boxes(bk).bg_gmm_pi) || isempty(up_boxes(bk).fg_gmm_pi)
            up_boxes(bk).fg_gmm_mu  = zeros(3,CP.NUM_GMM_CLUSTERS) + nan;
            up_boxes(bk).fg_gmm_cov = zeros(3,CP.NUM_GMM_CLUSTERS) + nan;
            up_boxes(bk).fg_gmm_pi  = zeros(CP.NUM_GMM_CLUSTERS,1) + nan;
            up_boxes(bk).bg_gmm_mu  = zeros(3,CP.NUM_GMM_CLUSTERS) + nan;
            up_boxes(bk).bg_gmm_cov = zeros(3,CP.NUM_GMM_CLUSTERS) + nan;
            up_boxes(bk).bg_gmm_pi  = zeros(CP.NUM_GMM_CLUSTERS,1) + nan;
        end     
      end
      [box_fg_2, box_bg_2] = eval_gmm_bboxes_mex(i1_bflt, up_boxes);
  
      % make selection of historic vs updated model
      for bk = 1:n_boxes;
        fg1 = box_fg_1(bk).prob > FGTHRESH;
        fg2 = box_fg_2(bk).prob > FGTHRESH;
        
        nfg1 = sum(fg1(:));
        nfg2 = sum(fg2(:));
        
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
    if exist('pi_bg','var')
        bad_gmm_bg = isnan(sum(pi_bg, 1));
    else
        bad_gmm_bg = [];
    end
    if exist('pi_fg','var')
        bad_gmm_fg = isnan(sum(pi_fg, 1));
    else
        bad_gmm_fg = [];
    end
    bad_gmm = bad_gmm_fg | bad_gmm_bg;
    all_boxes(bad_gmm) = [];
    
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
  end
end
end
