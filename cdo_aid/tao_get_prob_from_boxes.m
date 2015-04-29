function [boxes, prob_box_fg, prob_box_bg, count_box_fg, count_box_bg, weights_box] = tao_get_prob_from_boxes( ...
  occb_mask, past, i1_bflt, boxes, opts)

[rows, cols] = size(occb_mask);
imsize = [rows, cols];
% v2struct(opts);

if ~isempty(past.layers) && opts.CAUSAL && opts.BOXHELP && ~isempty(boxes);
  %----------------------------------------------------------
  % warp forward boxes
  % out of date: b.fg_prob{,_colour,_shape}
  %----------------------------------------------------------
  % [boxes_new, valid] = warp_bboxes_forward(boxes, past.layers, uvb_rev_bflt);
  [boxes_new, ~] = warp_bboxes_forward(boxes, past.layers, past.uvf_cbf_bflt);
  n_boxes_new = size(boxes_new, 1);
  
  %----------------------------------------------------------
  % compute new mask at current time
  %----------------------------------------------------------
  % evaluate colour
  [box_fg, box_bg] = eval_gmm_bboxes_mex(i1_bflt, boxes_new);
  for bk = 1:n_boxes_new;
    boxes_new(bk).fg_prob_colour = box_fg(bk).prob;
    boxes_new(bk).bg_prob_colour = box_bg(bk).prob;
  end
  
  % evaluate shape
  [shape_fg, shape_bg, shape_conf, shape_sigma] ...
    = obtain_bbox_shapes(past.t0.layers, boxes_new);
  for bk = 1:n_boxes_new;
    boxes_new(bk).sigma_shape   = shape_sigma(bk);
    boxes_new(bk).conf_shape    = shape_conf(bk).conf;
    boxes_new(bk).fg_prob_shape = shape_fg(bk).prob;
    boxes_new(bk).bg_prob_shape = shape_bg(bk).prob;
  end
    
  % combine colour and shape and flow
  for bk = 1:n_boxes_new;
    bb = boxes_new(bk); 
    sz_conf_shape = size(bb.conf_shape);
    
    sz_fg_shape = size(bb.fg_prob_shape);
    sz_fg_colour = size(bb.fg_prob_colour);
    if all(sz_conf_shape == sz_fg_shape) && all(sz_conf_shape == sz_fg_colour);
      % use colour and shape
      boxes_new(bk).fg_prob = bb.conf_shape .* bb.fg_prob_shape ...
        + (1.0 - bb.conf_shape) .* bb.fg_prob_colour;
      
      % % colour only
      % boxes_new(bk).fg_prob = bb.fg_prob_colour * bb.conf_colour;
      boxes_new(bk).conf = bb.conf_colour;
    else
      fprintf('%s: fg problem\n\n', mfilename);
    end
      
    sz_bg_shape = size(bb.bg_prob_shape);
    sz_bg_colour = size(bb.bg_prob_colour);
    if all(sz_conf_shape == sz_bg_shape) && all(sz_conf_shape == sz_bg_colour);
      % colour and shape
      boxes_new(bk).bg_prob = bb.conf_shape .* bb.bg_prob_shape ...
        + (1.0 - bb.conf_shape) .* bb.bg_prob_colour;
      
      % % color only
      % boxes_new(bk).bg_prob = bb.bg_prob_colour * bb.conf_colour;
      boxes_new(bk).conf = bb.conf_colour;
    else
      fprintf('%s: bg problem\n\n', mfilename);
    end
  end
    
  % % % % combine colour and shape and flow
  % % % for bk = 1:n_boxes_new;
  % % %   bb = boxes_new(bk);
  % % %
  % % %   sz_conf_shape = size(bb.conf_shape);
  % % %
  % % %   sz_fg_shape = size(bb.fg_prob_shape);
  % % %   sz_fg_colour = size(bb.fg_prob_colour);
  % % %   if all(sz_conf_shape == sz_fg_shape) && all(sz_conf_shape == sz_fg_colour);
  % % %     boxes_new(bk).fg_prob = bb.conf_shape .* bb.fg_prob_shape ...
  % % %       + (1.0 - bb.conf_shape) .* bb.fg_prob_colour;
  % % %   else
  % % %     fprintf('\nfg problem\n\n');
  % % %   end
  % % %
  % % %   sz_bg_shape = size(bb.fg_prob_shape);
  % % %   sz_bg_colour = size(bb.fg_prob_colour);
  % % %   if all(sz_conf_shape == sz_bg_shape) && all(sz_conf_shape == sz_bg_colour);
  % % %     boxes_new(bk).bg_prob = bb.conf_shape .* bb.bg_prob_shape ...
  % % %       + (1.0 - bb.conf_shape) .* bb.bg_prob_colour;
  % % %   else
  % % %     fprintf('\nbg problem\n\n');
  % % %   end
  % % % end
    
  %----------------------------------------------------------
  % use boxes to affect optimization
  %
  % TODO: use box_conf to affect weights computation (more
  %   weight to intensity or more weight to flow)
  %----------------------------------------------------------
  bb_conf = cat(1, boxes_new.conf_colour);
  bb_conf_inter = cat(1, boxes_new.conf);
  lay_occd = cat(1, boxes_new.lay_occd);
  [prob_box_fg, count_box_fg] = bboxes_to_prob_image( ...
    imsize, boxes_new, 1, bb_conf, bb_conf_inter);
  % [prob_box_local_bg, count_box_local_bg] = bboxes_to_prob_image( ...
  %   imsize, boxes_new, 0, bb_conf, bb_conf_inter);
  [prob_box_bg, count_box_bg] = bboxes_to_prob_image(imsize, ...
    boxes_new, 0, bb_conf .* max(0.0, 1.0 - lay_occd), bb_conf_inter);
  
  if any(isnan(prob_box_bg) | isnan(prob_box_fg));
    fprintf('\nbad prob_box_images\n\n'); 
    prob_box_bg(isnan(prob_box_bg))=0.5;
    prob_box_fg(isnan(prob_box_fg))=0.5;
  end
  
  [w_bg, cnt_bg] = bboxes_to_pixelwise_weights(imsize, boxes_new, 0);
  [w_fg, cnt_fg] = bboxes_to_pixelwise_weights(imsize, boxes_new, 1);
  % TODO: a better function than this?
  weights_box = w_bg + w_fg - 1.0 * ((cnt_bg + cnt_fg) > 0);
  
  if any(isnan(weights_box(:)));
    fprintf('\nbad weights_bg or weights_fg\n\n'); 
  end
  
  boxes = boxes_new;
  
  % check updated boxes, make sure still moving forward
  if opts.VIS < 200;
    box_centers = bbox_to_box_center(boxes);
    [bimg, ~] = draw_boxes(past.t0.layers, box_centers, opts.BOX_RAD);
    fig(55); imagesc(bimg); drawnow; title(1);
    fig(56); imagesc(min(weights_box, [], 3)); drawnow; title(1);
  end
  
  % TODO: check/debug the boxes again
  if (opts.VIS < 160);
    fprintf('displaying boxes\n');
    startDisplayImage(diagboxes(boxes));
  end
else
  boxes = [];
  prob_box_fg = zeros(imsize);
  prob_box_bg = zeros(imsize);
  count_box_fg = zeros(imsize);
  count_box_bg = zeros(imsize);
  weights_box = zeros(imsize);
end
end
