function [out_img, box_img] = draw_boxes(img, boxes, sizes, conf)
out_img = img;
box_img = zeros(size(img));

if isempty(boxes);
  return;
end

[rows, cols, chan] = size(img);
imsize = [rows, cols];
[N, ~] = size(boxes);

if ~exist('sizes', 'var'); sizes = 5; end;
if length(sizes) == 1; sizes = sizes * ones(N, 1); end

box_img = zeros(imsize);
box_draw_img = zeros(imsize);
boxes = round(boxes);
if ~exist('conf', 'var'); conf = ones(N, 1); end;
for k = 1:N;
  binds = put_box(imsize, boxes(k, :), sizes(k));

  box_img(binds) = k;
  box_img(boxes(k,:)) = k;

  box_draw_img(binds) = conf(k);
  box_draw_img(boxes(k, :)) = conf(k);
end

if chan == 3;
  if isa(img, 'uint8');
    out_img = 0.6 * uint8(255.0 * sc(box_draw_img, 'jet')) + 0.4 * img;
  elseif isa(img, 'double');
    out_img = 0.6 * sc(box_draw_img, 'jet') + 0.4 * img;
  else
    fprintf('broke\n'); 
  end
elseif chan == 1;
  out_img = 0.6 * box_draw_img * max(img(:)) + 0.4 * img;
else
  fprintf('broke\n');
end
end

% returns the image indices of the pixels to be drawn on
function draw_box_inds = put_box(imsize, box, sz)
rows = imsize(1);
cols = imsize(2);

y = box(1);
x = box(2);

ymin = max(1,    y - sz);
ymax = min(rows, y + sz);
xmin = max(1,    x - sz);
xmax = min(cols, x + sz);

ys = (ymin:ymax)';
xs = (xmin:xmax)';

nys = length(ys);
nxs = length(xs);

sy = ones(nys-1, 1);
sx = ones(nxs-1, 1);
top   = [ys(1) * sx  , xs(2:end)   ];
bot   = [ys(end) * sx, xs(1:end-1) ];
left  = [ys(1:end-1) , xs(1) * sy  ];
right = [ys(2:end)   , xs(end) * sy];

draw_box_yx = cat(1, top, right, bot, left);
draw_box_inds = sub2ind(imsize, draw_box_yx(:, 1), draw_box_yx(:, 2));
end
