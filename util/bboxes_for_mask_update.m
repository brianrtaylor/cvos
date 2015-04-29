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
function [new_boxes, boxes, valid_boxes] = bboxes_for_mask_update( ...
  mask_in, Dx, Dy, r, boxes, box_limit, boxes_per_pixel_limit)

BOX_MEMORY_FACTOR = 0.2;
VIS = 999999999;

if ~exist('box_limit', 'var') || (box_limit < 0); box_limit = inf; end;
if ~exist('boxes_per_pixel_limit', 'var'); boxes_per_pixel_limit = 3; end;
[rows, cols, ~] = size(mask_in);
imsize = [rows, cols];
dsk = strel('disk', 5);

% mask = double(mask_in(:) > 0.5);
mask = round(mask_in);
mask_vals = unique(mask_in);
n_mask_vals = length(mask_vals);

if isempty(boxes); n_boxes = 0; else; n_boxes = size(boxes, 1); end

% now per box item
% [box_height, box_width, ~] = size(box_fg);
diam = 2 * r + 1;
% npx = diam * diam;
    
yy = repmat((-r:r)', 1, diam);
xx = repmat((-r:r) , diam, 1);
box_dist = sqrt(yy .* yy + xx .* xx);
inv_box_dist = 1./ max(1.0, box_dist);

%--------------------------------------------------------------------
% build edge image to cover
%--------------------------------------------------------------------
% weights edge weight
dx_l = Dx * mask(:);
dy_l = Dy * mask(:);
dx_l = reshape(dx_l, imsize);
dy_l = reshape(dy_l, imsize);
d_l  = dx_l + dy_l;
ad_l = abs(d_l);

look_img = abs(d_l) > 0;
look = look_img;

%--------------------------------------------------------------------
% go through existing boxes and choose to keep / toss
%--------------------------------------------------------------------
% mark first, delete in fell swoop after
% p_all = find(look_img > 0);
% check_p = look(p_all);
counts = zeros(imsize);

% methods to balance boxes:
% 1. have a list per pixel. when we count too much, remove 1 from the pixels
%    and also remove the that box from the list of blocks per pixel (way slow)
%    stop when no pixel has more than 5 or so on it
%
% 2. sort the boxes in advance by their confidence. once certain pixels hit 
%    5 or so on it, don't add any more boxes that might touch this pixel. 
%    Just ignore them and drop from the box set. This is probably the best
%    option since it's of reasonable speed. now try implementing it
%    tomorrow. a. do that sorting by confidence and b. go

% sortboxes on confidence
if ~isempty(boxes);
  [conf_bb_sorted, boxes_ind] = sort(cat(1, boxes.conf), 1, 'descend');
  boxes = boxes(boxes_ind);
  
  if VIS < 100;
    if ~isempty(conf_bb_sorted);
      fprintf('\n--- boxes conf ---\n');
      stat(conf_bb_sorted);
      fprintf('\n\n');
    end
  end
end

valid_boxes = true(n_boxes, 1);
for k = 1:n_boxes;
  bb = boxes(k);
  
  y = round(boxes(k).y);
  x = round(boxes(k).x);
  r = round(boxes(k).r);
  diam = 2 * r + 1;
  z = zeros(diam, diam);
  
  BOXJITTER = 1;
  if BOXJITTER;
    % modify y,x based on the edge image
    ymin = max(1,    y - r);
    ymax = min(rows, y + r);
    xmin = max(1,    x - r);
    xmax = min(cols, x + r);
    
    ys = ymin:ymax;
    xs = xmin:xmax;
    
    by = y - ymin + 1;
    bx = x - xmin + 1;
    
    % TODO: can use boxes(k).{x,y} here instead of x,y if you'd like
    bedge = ad_l(ys,xs) > 0;
    [bedgey, bedgex] = find(bedge);
    if ~isempty(bedgey) && ~isempty(bedgex);
      bdist = (bedgey - boxes(k).y) .^ 2 + (bedgex - boxes(k).x) .^ 2;
      [~, eidx] = min(bdist);
      y = ymin + bedgey(eidx) - 1;
      x = xmin + bedgex(eidx) - 1;
%     else
%       % stop early check
%       boxes(k).conf = boxes(k).conf  - BOX_MEMORY_FACTOR;
%       % valid_boxes(k) = false; continue;
    end
  end
  
  % image
  ymin = max(1,    y - r);
  ymax = min(rows, y + r);
  xmin = max(1,    x - r);
  xmax = min(cols, x + r); 
  
  ys = ymin:ymax;
  xs = xmin:xmax;

  nys = length(ys);
  nxs = length(xs);
  npx = nys*nxs;
  
  % box
  bymin = ymin - y + r + 1;
  bymax = ymax - y + r + 1;
  bxmin = xmin - x + r + 1;
  bxmax = xmax - x + r + 1;
  
  bys = bymin:bymax;
  bxs = bxmin:bxmax;

%   % fg/bg prob
%   fg_prob = box_fg(:,:,k);
%   bg_prob = box_bg(:,:,k);
  
  %------------------------------------------------------------------
  % changes to the box
  %
  % TODO: change this portion to preserve the box as it was before if the
  % layer has changed to 0, cuz then maybe we lost the object, but it's
  % still there. sooo.. check that
  %------------------------------------------------------------------
  % udpate age
  boxes(k).age = boxes(k).age + 1;
  
  % update lay_occd, lay_occr
  lay_patch = z;
  lay_patch(bys,bxs) = mask(ys, xs);
  
  patch_hist_fg = zeros(n_mask_vals, 1);
  patch_hist_bg = zeros(n_mask_vals, 1);
  for vk = 1:n_mask_vals;
    patch_hist_fg(vk) = sum(vec(bb.fg_prob(lay_patch == mask_vals(vk))));
    patch_hist_bg(vk) = sum(vec(bb.bg_prob(lay_patch == mask_vals(vk))));
  end
  
  [~, fg_vote] = max(patch_hist_fg);
  fg_choice = mask_vals(fg_vote);
  if fg_choice > 0;
    boxes(k).lay_occr = fg_choice;
    boxes(k).lay_occr_conf = patch_hist_fg(fg_vote) / sum(patch_hist_fg);

    [~, bg_vote] = max(patch_hist_bg);
    boxes(k).lay_occd = mask_vals(bg_vote);
    boxes(k).lay_occd_conf = patch_hist_bg(bg_vote) / sum(patch_hist_bg);

    % update foreground background masks
    boxes(k).fg_mask = double(abs(lay_patch - boxes(k).lay_occr) < 0.5);
    boxes(k).bg_mask = double(abs(lay_patch - boxes(k).lay_occd) < 0.5);
    boxes(k).gmm_learn_colour_mask = double(1 ...
      - imdilate(edge(boxes(k).fg_mask), dsk));
  else
    boxes(k).conf = boxes(k).conf - BOX_MEMORY_FACTOR;
    % boxes(k).conf_colour = boxes(k).conf_colour * 0.75;
  end
  
  % update patches
  % boxes(k).patches = cat(3, boxes(k).patches, i1(ys, xs, :));

  %------------------------------------------------------------------
  % keep or toss box criteria
  %------------------------------------------------------------------  
  % box outside of image limits
  if (y < 1) || (y > rows) || (x < 1) || (x > cols);
    valid_boxes(k) = false; continue;
  end

  % % box shrinks too much (image edge)
  % if ((nys * nxs) < (0.5 * npx));
  %   valid_boxes(k) = false; continue;
  % end

  % if fg or bg < %10 of the box, lower confidence
  % or if the mask to learn gmm from / do movement is too small
  msk = boxes(k).fg_mask > 0.5;
  msk_learn = msk .* boxes(k).fg_prob .* boxes(k).gmm_learn_colour_mask;
  nmsk = sum(msk(:));
  nmsk_learn = sum(msk_learn(:));
  TOOFEW = 0.1 * npx;
  TOOMUCH = 0.9 * npx;
  if (nmsk < TOOFEW) || (nmsk > TOOMUCH) || (nmsk_learn < TOOFEW) ...
    || (nmsk_learn > TOOMUCH);
    boxes(k).conf = boxes(k).conf - BOX_MEMORY_FACTOR;
    % boxes(k).conf_colour = boxes(k).conf_colour - BOX_MEMORY_FACTOR;
    % valid_boxes(k) = false; continue;
  end
  
  % if new lay_occr is background layer, or occlusion constraint isn't
  % maintained, then toss
  if (boxes(k).lay_occr == 0) || (boxes(k).lay_occr <= boxes(k).lay_occd);
    valid_boxes(k) = false; continue;
  end
  
  % if confidences in colour and shape fall too low
  if boxes(k).conf_colour < 0.25;
    valid_boxes(k) = false; continue;
  end
 
  % general confidence too low
  if boxes(k).conf < 0.10;
    valid_boxes(k) = false; continue;
  end

  % too many boxes in one spot
  if any(vec(counts(ys, xs))) >= boxes_per_pixel_limit;
    valid_boxes(k) = false; continue;
  end

  %------------------------------------------------------------------
  % if kept, edit the look image for box addition step
  %------------------------------------------------------------------
  look(ys, xs) = 0;

  % for seeing if too many boxes
  counts(ys, xs) = counts(ys, xs) + 1;
end

if ~isempty(valid_boxes);
  boxes = boxes(valid_boxes);
  n_boxes = sum(valid_boxes);
end

if n_boxes > 0;
  if VIS < 100;
    stat(cat(1, boxes.conf));
    fprintf('\n--------------------\n');
  end
end

%--------------------------------------------------------------------
% add new boxes
%--------------------------------------------------------------------
p_all = find(look > 0);
% [py, px] = ind2sub(imsize, p_all);

p_out       = zeros(1, 0);
pts_to_add  = zeros(1, 0);
check_p     = look(p_all);

new_boxes = struct;
new_boxes.y = [];
new_boxes.x = [];
new_boxes.r = [];
new_boxes.d = [];
new_boxes.invd = [];
new_boxes.lay_occd = [];
new_boxes.lay_occr = [];
new_boxes.lay_occd_conf = [];
new_boxes.lay_occr_conf = [];
new_boxes.fg_mask = []; % from layers (might be old)
new_boxes.bg_mask = [];
new_boxes.fg_prob = []; % from colour + shape classification
new_boxes.bg_prob = [];
new_boxes.fg_prob_colour = [];
new_boxes.bg_prob_colour = [];
new_boxes.fg_prob_shape = [];
new_boxes.bg_prob_shape = [];
new_boxes.gmm_change = [];
new_boxes.gmm_learn_colour_mask = [];
new_boxes.fg_gmm_mu = [];
new_boxes.fg_gmm_cov = [];
new_boxes.fg_gmm_pi = [];
new_boxes.bg_gmm_mu = [];
new_boxes.bg_gmm_cov = [];
new_boxes.bg_gmm_pi = [];
new_boxes.conf = [];
new_boxes.conf_colour = [];
new_boxes.conf_colour_denom = [];
new_boxes.conf_colour_h = [];
new_boxes.conf_colour_denom_h = [];
new_boxes.conf_mask = [];
new_boxes.conf_shape = [];
new_boxes.sigma_shape = [];
new_boxes.age = [];
new_boxes.patches = [];
new_boxes.u_bg = [];
new_boxes.v_bg = [];
new_boxes.u_fg = [];
new_boxes.v_fg = [];
new_boxes.conf_uv_bg = []; % some variance
new_boxes.conf_uv_fg = []; % some variance
new_boxes.bg_prob_uv = []; % some variance
new_boxes.fg_prob_uv = []; % some variance

k = 0;
nboxes = 0;
MAX_PIXELS_TO_IGNORE = 2;
while sum(check_p) > MAX_PIXELS_TO_IGNORE;
  % testing, enough boxes to deal with for now
  if size(p_out, 1) >= (box_limit - n_boxes);
    break;
  end

  if ~isempty(pts_to_add);
    % extract a point and add to output list
    p = pts_to_add(1);
    pts_to_add = pts_to_add(2:end);

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
      % fprintf('skip (%d, %d)\n', py, px);
      continue;
    end

    %-------------------------------------------------
    % get more box info
    %-------------------------------------------------
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
    
    assert(occluder > occluded, sprintf( ...
      '%s l.350: occluder/occluded selection failing', mfilename));
    
    %-------------------------------------------------
    % make box and add it
    %-------------------------------------------------   
    p_out = [p_out; p];
   
    nboxes = nboxes + 1;
    bb_new = struct;
    bb_new.y = double(py);
    bb_new.x = double(px);
    bb_new.r = double(r);
    bb_new.d = box_dist;
    bb_new.invd = inv_box_dist;
    bb_new.lay_occd = double(occluded);
    bb_new.lay_occr = double(occluder);
    bb_new.lay_occd_conf = 1.0;
    bb_new.lay_occr_conf = 1.0;
    bb_new.fg_mask = double(abs(mask_in(ys,xs) - occluder) < 0.5); 
    bb_new.bg_mask = double(abs(mask_in(ys,xs) - occluded) < 0.5); 
    bb_new.fg_prob = bb_new.fg_mask;
    bb_new.bg_prob = bb_new.bg_mask;
    bb_new.fg_prob_colour = double([]);
    bb_new.bg_prob_colour = double([]);
    bb_new.fg_prob_shape = double([]);
    bb_new.bg_prob_shape = double([]);
    
    bb_new.gmm_change = double([]);
    bb_new.gmm_learn_colour_mask = double(1 ...
      - imdilate(edge(bb_new.fg_mask), dsk));
    bb_new.fg_gmm_mu = double([]);
    bb_new.fg_gmm_cov = double([]);
    bb_new.fg_gmm_pi = double([]);
    bb_new.bg_gmm_mu = double([]);
    bb_new.bg_gmm_cov = double([]);
    bb_new.bg_gmm_pi = double([]);
    
    bb_new.conf = double(1.0);
    bb_new.conf_colour = double([]);
    bb_new.conf_colour_denom = double([]);
    bb_new.conf_colour_h = double([]);
    bb_new.conf_colour_denom_h = double([]);
    bb_new.conf_mask = double([]);
    bb_new.conf_shape = double([]);
    bb_new.sigma_shape = double([]);
    
    bb_new.age = double(0);
    bb_new.patches = double([]);
    % % if exist('i1', 'var') && ~isempty(i1) && 0; % don't do for now
    % %   bb_new.patches = double(i1(ys,xs,:));
    % % end

    bb_new.u_bg = [];
    bb_new.v_bg = [];
    bb_new.u_fg = [];
    bb_new.v_fg = [];
    bb_new.conf_uv_bg = []; % some variance
    bb_new.conf_uv_fg = []; % some variance
    bb_new.bg_prob_uv = [];
    bb_new.fg_prob_uv = [];
    
    new_boxes(nboxes, 1) = bb_new;

    
    % find new points and add them to make boxes for
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
% [py, px] = ind2sub(imsize, p_out);
% boxes_new = [py, px];
% n_new = size(p_out, 1);

%--------------------------------------------------------------------
% output
%--------------------------------------------------------------------
% boxes       = [boxes; boxes_new];
% layer_occd  = [lay_occd; new_layer_occd];
% layer_occr  = [lay_occr; new_layer_occr];
% box_ages    = [(box_ages + 1); zeros(n_new, 1)];
% box_conf    = [box_conf; ones(n_new, 1)];
% 
% assert(size(boxes, 1) == size(layer_occd, 1) ...
%   && size(boxes, 1) == size(layer_occr, 1), 'boxes size ~= layer labels size');
if nboxes == 0; new_boxes = []; end
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
