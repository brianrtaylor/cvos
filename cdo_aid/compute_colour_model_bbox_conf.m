function [conf, denom] = compute_colour_model_bbox_conf(boxes, layers)

n_boxes = length(boxes);
[rows, cols, ~] = size(layers);
imsize = [rows, cols];

conf = zeros(n_boxes, 1);
denom = zeros(n_boxes, 1);

for k = 1:n_boxes;
  bb = boxes(k);
  y = bb.y;
  x = bb.x;
  r = bb.r;
  
  % image
  ymin = max(1,    y - r);
  ymax = min(rows, y + r);
  xmin = max(1,    x - r);
  xmax = min(cols, x + r);
  
  ys = ymin:ymax;
  xs = xmin:xmax;

  % box
  bymin = ymin - y + r + 1;
  bymax = ymax - y + r + 1;
  bxmin = xmin - x + r + 1;
  bxmax = xmax - x + r + 1;
  
  bys = bymin:bymax;
  bxs = bxmin:bxmax;
  
  % lay_patch = layers(ys,xs);
  % fg_mask = abs(lay_occd(k) - lay_patch)  > abs(lay_occr(k) - lay_patch);   
  fg_edge = double(edge(bb.fg_mask));
  [ii, jj] = find(fg_edge);
  a = msfm(1 - fg_edge, [ii, jj]', true, true);
  weight_dist = exp(- (a.*a) ./ (r*r));
  % weight_dist = ones(size(lay_patch));
  box_prob = bb.fg_prob;

  try
    weighted_fg_wrong = (abs(bb.fg_mask - box_prob) .* weight_dist);
  catch e
    e
  end
  
  if bb.gmm_change == true;
    if bb.age == 0; % isempty(bb.conf_colour);
      denom(k) = 1e-8 + sum(weight_dist(:));
      conf(k) = 1 - (sum(weighted_fg_wrong(:)) / denom(k));
    else
      denom(k) = bb.conf_colour_denom_h + sum(weight_dist(:));
      old_weighted_fg_wrong = (1 - bb.conf_colour_h) * bb.conf_colour_denom_h;
      conf(k) = 1 - (old_weighted_fg_wrong + sum(weighted_fg_wrong(:))) / denom(k);
      % denom(k) = bb.conf_colour_denom_h + sum(weight_dist(:));
      % old_weighted_fg_wrong = (1 - bb.conf_colour_h) * bb.conf_colour_denom_h;
      % conf(k) = 1 - (old_weighted_fg_wrong + sum(weighted_fg_wrong(:))) / denom(k);
    end
  else
    conf(k) = bb.conf_colour;
  end
  
  if conf(k) == 0;
    fprintf('\nzero conf colour\n\n');
  end
end
end
