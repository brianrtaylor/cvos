function [weights, counts, weightsx, weightsy, countsx, countsy] = boxes_to_pixelwise_weights( ...
  imsize, boxes, boxes_prob, box_conf)

n_boxes = size(boxes, 1);
boxes = round(boxes);
if ~exist('box_conf', 'var'); box_conf = ones(n_boxes, 1); end;
weightsx = zeros(imsize);
weightsy = zeros(imsize);
countsx = zeros(imsize);
countsy = zeros(imsize);

rows = imsize(1);
cols = imsize(2);

[height, width, ~] = size(boxes_prob);
r = (height - 1) / 2;

for k = 1:n_boxes;
  y = boxes(k, 1);
  x = boxes(k, 2);
  
  % image domain 
  ymin = max(1,    y - r);
  ymax = min(rows, y + r);
  xmin = max(1,    x - r);
  xmax = min(cols, x + r);
  
  ys = ymin:ymax;
  xs = xmin:xmax;

  % box domain
  bymin = ymin - y + r + 1;
  bymax = ymax - y + r + 1;
  bxmin = xmin - x + r + 1;
  bxmax = xmax - x + r + 1;

  bheight = bymax - bymin + 1;
  bwidth  = bxmax - bxmin + 1;
  
  bys = bymin:bymax;
  bxs = bxmin:bxmax;

  box_prob = boxes_prob(bys,bxs,k);

  pdiffx = abs(box_prob(:,1:(end-1)) - box_prob(:, 2:end));
  pdiffx_box = reshape(pdiffx, [bheight, bwidth-1]);
  pdiffy = abs(box_prob(1:(end-1),:) - box_prob(2:end, :));
  pdiffy_box = reshape(pdiffy, [bheight-1, bwidth]);
  wx = 1 - pdiffx_box;
  wy = 1 - pdiffy_box;
  
  if any(vec(isnan(wx))) || any(vec(isnan(wy)));
    fprintf('NaNs inside of a box probability\n');
    wx(isnan(wx))=0;
    wy(isnan(wy))=0;
  end

  weightsx(ys, xs(1:end-1)) = weightsx(ys, xs(1:end-1)) + box_conf(k) * wx;
  weightsy(ys(1:end-1), xs) = weightsy(ys(1:end-1), xs) + box_conf(k) * wy;
  
  countsx(ys, xs(1:end-1)) = countsx(ys, xs(1:end-1)) + 1;
  countsy(ys(1:end-1), xs) = countsy(ys(1:end-1), xs) + 1;
  
end

weightsx = weightsx ./ max(1.0, countsx);
weightsy = weightsy ./ max(1.0, countsy);
weights  = [weightsx(:); weightsy(:)];
counts   = [countsx(:); countsy(:)];
end
