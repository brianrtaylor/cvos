function prob = boxes_to_prob_image(imsize, boxes, box_fg, box_conf)

n_boxes = size(boxes, 1);
if ~exist('box_conf', 'var'); box_conf = ones(n_boxes, 1); end;
% [rows, cols, chan] = size(i1);
prob = zeros(imsize);
count = ones(imsize);
rows = imsize(1);
cols = imsize(2);

boxes = round(boxes);
[height, width, ~] = size(box_fg);
r = (height - 1) / 2;

for k = 1:n_boxes;
  y = boxes(k, 1);
  x = boxes(k, 2);
  
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
  
  % combine via averaging
  prob(ys, xs) = prob(ys, xs) + box_fg(bys,bxs,k) .* box_conf(k);
  
  count(ys, xs) = count(ys, xs) + 1;
end

prob = prob ./ count;
end
