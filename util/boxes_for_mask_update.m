% @return: boxes: boxes are defined by a N x 3 array of centers (y,x), and
% a radius
%
% @note: how to treat boxes
% * decay over time, 0.2 or something?
% * remove if only 5% fg or only 5% bg
% * remove if the box shrinks (below an X%?)
%
% @note: things boxes need
% * center (y,x)
% * radii - or size of some sort
% * layer_occd
% * layer_occr
% * mask (perhaps)
% * colour model (box_gmm_bg, box_gmm_fg)
% * colour confidence
% * shape model 
function [boxes, layer_occr, layer_occd, box_bg, box_fg, box_ages, box_conf] = boxes_for_mask_update( ...
  mask_in, Dx, Dy, minspace, r, boxes, box_bg, box_fg, ...
  lay_occr, lay_occd, box_ages, box_conf, box_limit)

if ~exist('box_limit', 'var') || (box_limit < 0); box_limit = inf; end;
if ~isempty(boxes); boxes = round(boxes); end;
% szs = r - 1;
[rows, cols, ~] = size(mask_in);
imsize = [rows, cols];
mask = double(mask_in(:) > 0.5);
if isempty(boxes); n_boxes = 0; else; n_boxes = size(boxes, 1); end
[box_height, box_width, ~] = size(box_fg);
diam = 2 * r + 1;
npx = diam * diam;

%--------------------------------------------------------------------
% build edge image to cover
%--------------------------------------------------------------------
% weights edge weight
dx_l = Dx * mask;
dy_l = Dy * mask;
dx_l = reshape(dx_l, imsize);
dy_l = reshape(dy_l, imsize);
d_l  = dx_l + dy_l;

look_img = abs(d_l) > 0;
look = look_img;

%--------------------------------------------------------------------
% go through existing boxes and choose to keep / toss
%--------------------------------------------------------------------
% mark first, delete in fell swoop after
valid_boxes = true(n_boxes, 1);
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

  nys = length(ys);
  nxs = length(xs);
  
%   % box
%   bymin = ymin - y + r + 1;
%   bymax = ymax - y + r + 1;
%   bxmin = xmin - x + r + 1;
%   bxmax = xmax - x + r + 1;
  
%   bys = bymin:bymax;
%   bxs = bxmin:bxmax;

  % fg/bg prob
  fg_prob = box_fg(:,:,k);
  bg_prob = box_bg(:,:,k);

  %------------------------------------------------------------------
  % keep or toss box criteria
  %------------------------------------------------------------------
  % box outside of image limits
  if (y < 1) || (y > rows) || (x < 1) || (x > cols);
    valid_boxes(k) = false; continue;
  end

  % box shrinks too much (image edge)
  if ((nys * nxs) < (0.5 * npx));
    valid_boxes(k) = false; continue;
  end

  % if fg or bg < 10% of the box, lower confidence
  mask_fg = fg_prob > bg_prob;
  if (sum(mask_fg(:)) < (0.2 * npx)) || (sum(mask_fg(:)) > (0.8 * npx));
    box_conf(k) = box_conf(k) * 0.5;
    % valid_boxes(k) = false; continue;
  end
  
  % if confidences in colour and shape fall too low

  %------------------------------------------------------------------
  % if kept, edit the look image for box addition step
  %------------------------------------------------------------------
  look(ys, xs) = 0;
end

if ~isempty(valid_boxes);
  n_boxes = sum(valid_boxes);
  boxes = boxes(valid_boxes, :);
  box_fg = box_fg(:, :, valid_boxes);
  box_bg = box_bg(:, :, valid_boxes);
  lay_occd = lay_occd(valid_boxes);
  lay_occr = lay_occr(valid_boxes);
  box_ages = box_ages(valid_boxes);
end

%--------------------------------------------------------------------
% add new boxes
%--------------------------------------------------------------------
% while there are still not turned off edge pixels left
% TOOD: later change this to some percentage (instead of like any dot)
p_all = find(look_img > 0);
% [py, px] = ind2sub(imsize, p_all);

p_out       = zeros(1, 0);
pts_to_add  = zeros(1, 0);
check_p     = look(p_all);

new_layer_occd    = zeros(1, 0);
new_layer_occr    = zeros(1, 0);

k = 0;
while sum(check_p) > 0;
  % testing, enough boxes to deal with for now
  if size(p_out, 1) >= (box_limit - n_boxes);
    break;
  end

  if ~isempty(pts_to_add);
    % extract a point and add to output list
    p = pts_to_add(1);
    pts_to_add = pts_to_add(2:end);


    % make a box 

    % TODO: deal with ignoring boxes that go out of bounds, aka don't 
    % add them, but remove their pixels as desired. or maybe shift the 
    % box in bounds (maybe latter) --> for now, dropped
    [py, px] = ind2sub(imsize, p);
    ys = ((py - r):(py + r))';
    xs = ((px - r):(px + r))';
    if ((ys(1) < 1) || (ys(end) > rows) || (xs(1) < 1) || (xs(end) > cols));
      ys = ((max(1, py - r)):(min(rows, py + r)))';
      xs = ((max(1, px - r)):(min(cols, px + r)))';
      look(ys, xs) = 0;
      check_p = look(p_all);
%       fprintf('skip (%d, %d)\n', py, px);
      continue;
    end

    p_out = [p_out; p];
    edge_inds = put_box(imsize, [py, px], r);
    
    % get layer values for occluder, occluded
    if dx_l(p) > 0;
      occluded = mask_in(p);
      occluder = mask_in(p + rows);
    elseif dx_l(p) < 0;
      occluded = mask_in(p + rows);
      occluder = mask_in(p);
    elseif dy_l(p) > 0;
      occluded = mask_in(p);
      occluder = mask_in(p + 1);
    elseif dy_l(p) < 0;
      occluded = mask_in(p + 1);
      occluder = mask_in(p);
    end
    
    new_layer_occd = [new_layer_occd; occluded];
    new_layer_occr = [new_layer_occr; occluder];

    % find new points and add them
    p_new = edge_inds(look(edge_inds) > 0);
    % temp gross hack
    if ~isempty(p_new);
      p_new = p_new(1);
      % pts_to_add = cat(1, pts_to_add, p_new(:));
      pts_to_add = [pts_to_add; p_new(:)];
    end

    % draw on look
    look(ys, xs) = 0;
    check_p = look(p_all);
  else
    pts_to_add = p_all(find(check_p > 0, 1)); % take the first unset one
  end
  k = k + 1;

  DRAW = false;
  if ~isempty(p_out) && DRAW;
    [py, px] = ind2sub(imsize, p_out);
    p_to_draw = [py, px];
    [out_img, box_img] = draw_boxes(mask_in, p_to_draw, r);
    fig(55); clf; imagesc(out_img);
    drawnow; pause(0.1);
  end
end

% TODO: before we do this step, need to go from some p_all to some ordered
% list of items, ordered in the sense of going in a continuous line around
% the object. I think this means we're gonna do some kind of walk algorithm

% p = p_all(1:minspace:length(p_all));
[py, px] = ind2sub(imsize, p_out);
boxes_new = [py, px];
n_new = size(p_out, 1);

%--------------------------------------------------------------------
% output
%--------------------------------------------------------------------
boxes       = [boxes; boxes_new];
layer_occd  = [lay_occd; new_layer_occd];
layer_occr  = [lay_occr; new_layer_occr];
box_ages    = [(box_ages + 1); zeros(n_new, 1)];

% TODO: remove confidence over time?
% TODO: remove boxes of too low confidence
box_conf    = [box_conf; ones(n_new, 1)];

assert(size(boxes, 1) == size(layer_occd, 1) ...
  && size(boxes, 1) == size(layer_occr, 1), 'boxes size ~= layer labels size');
end

%--------------------------------------------------------------------
% helper functions
%--------------------------------------------------------------------
function edge_inds = put_box(imsize, box, sz)

y = box(1);
x = box(2);

ys = ((y - sz):(y + sz))';
xs = ((x - sz):(x + sz))';

n = 2 * sz + 1;
s = ones(n-1, 1);
top   = [ys(1) * s, xs(2:end)  ];
bot   = [ys(n) * s, xs(1:end-1)];
left  = [ys(1:end-1), xs(1) * s];
right = [ys(2:end)  , xs(n) * s];

draw_box_yx = cat(1, top, right, bot, left);
edge_inds = sub2ind(imsize, draw_box_yx(:, 1), draw_box_yx(:, 2));
end
