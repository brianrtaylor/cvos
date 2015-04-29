function [alpha_map_bg, alpha_map_fg] = boxes_to_alpha(imsize, layers, ...
  layer_occr, layer_occd, boxes, box_bg, box_fg)
% TODO: write a function to create an alpha map with the following primary function.
% essentially copy the other box iteration mfunctions you got


%         box_alpha = params.PROB_BOX_FG ...
%           * ((layer_occr(k) == past.t0.layers) .* prob_box_fg ...
%           +  (layer_occd(k) == past.t0.layers) .* prob_box_bg);

boxes = round(boxes);
n_boxes = size(boxes, 1);
alpha_map_fg = zeros(imsize);
alpha_map_bg = zeros(imsize);
counts_fg = ones(imsize);
counts_bg = ones(imsize);
rows = imsize(1);
cols = imsize(2);
diam = size(box_bg, 1);
r = (diam - 1) / 2;

for k = 1:n_boxes;
  y = boxes(k, 1);
  x = boxes(k, 2);
  
  ymin = max(1,    y - r);
  ymax = min(rows, y + r);
  xmin = max(1,    x - r);
  xmax = min(cols, x + r);
  
  ys = ymin:ymax;
  xs = xmin:xmax;
  
  lay_patch = layers(ys,xs);
  
  fg_mask = abs(layer_occd(k) - lay_patch)  > abs(layer_occr(k) - lay_patch); 
  alpha_add_fg = fg_mask .* box_fg(:,:,k);
  alpha_map_fg(ys, xs) = alpha_map_fg(ys, xs) + alpha_add_fg;
  alpha_add_bg = (1 - fg_mask) .* box_bg(:,:,k);
  alpha_map_bg(ys, xs) = alpha_map_bg(ys, xs) + alpha_add_bg;
  
  % TODO: later replace with conf
  counts_fg(ys, xs) = counts_fg(ys, xs) + 1;
  counts_bg(ys, xs) = counts_bg(ys, xs) + 1;
  
end

alpha_map_fg = alpha_map_fg ./ counts_fg;
alpha_map_bg = alpha_map_bg ./ counts_bg;
end