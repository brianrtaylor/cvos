function [shape_fg, shape_bg, shape_conf, sigma] = obtain_bbox_shapes( ...
  layers, boxes)

n_boxes = size(boxes, 1);
[rows, cols, ~] = size(layers);

shape_fg = struct('prob', []);
shape_bg = struct('prob', []);
shape_conf = struct('conf', []);
sigma = nan(n_boxes, 1);

% parameters
CONF_COLOUR_THRESH = 0.85;
SIGMA_LO = 2.0;

a_denom = (1 - CONF_COLOUR_THRESH) ^ 2;

for k = 1:n_boxes;
  bb = boxes(k);
  
  y = round(bb.y);
  x = round(bb.x);
  r = round(bb.r);
  diam = 2 * r + 1;
  z = zeros(diam, diam);
  
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
  
  out_of_bounds = 1 - z;
  z(bys, bxs) = 0;
 
  % 1. find new shape mask
  lay_patch = z;
  lay_patch(bys, bxs) = layers(ys, xs);
  shape_fg(k).prob = double(abs(lay_patch - bb.lay_occr) < 0.5);
  shape_bg(k).prob = double(abs(lay_patch - bb.lay_occd) < 0.5);
  shape_fg(k).prob(out_of_bounds) = 0;
  shape_bg(k).prob(out_of_bounds) = 0;
  
  % 1. compute shape sigma
  if bb.conf_colour < CONF_COLOUR_THRESH;
    sigma(k) = SIGMA_LO;
  else
    a = (diam - SIGMA_LO) / a_denom;
    sigma(k) = SIGMA_LO + a * (bb.conf_colour - CONF_COLOUR_THRESH)^2;
  end
  
  % 2. compute shape conf
  fg_edge = double(edge(shape_fg(k).prob));
  [ii, jj] = find(fg_edge);
  a = msfm(1.0 - fg_edge, [ii, jj]', true, true);
  shape_conf(k).conf = 1.0 - exp(- (a.*a) ./ (sigma(k) * sigma(k)));
end
end
