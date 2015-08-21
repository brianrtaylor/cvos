%-----------------------------------------------------------------------------
% bboxes_for_mask_update
%
% updates the set of local shape classifiers by removing those that poorly 
% approximate the object boundary and adding new ones where the boundary is
% not yet covered by any. A box represents a local shape classifier
% 
% @return: new_boxes: struct array of the boxes created 
% @return: boxes: struct array of the updated set of boxes
% @return: valid_boxes: indicates which boxes are still good and which should
%   be dropped (based on criteria including how much foreground is in the 
%   classifier window. Search "keep or toss box criteria" in this file
%
% @param: mask_in (MxN): input object mask
% @param: Dx, Dy: difference operators in x, y directions
% @param: r: radius for new boxes
% @param: boxes: struct array of current boxes
% @param: box_limit: limit on total number of boxes (default: infinite)
% @param: boxes_per_pixel_limit: max number of boxes covering each pixel
%   (if too many cover a pixel, the lowest model confidence are dropped first)
%   (default: 3)
%
% @note: box structure description
% * center (y,x)
% * radii - or size of some sort
% * layer_occd
% * layer_occr
% * mask (perhaps)
% * colour model (box_gmm_bg, box_gmm_fg)
% * colour confidence
% * shape model
%-----------------------------------------------------------------------------
function [new_boxes, boxes, valid_boxes] = bboxes_for_mask_update( ...
  mask_in, Dx, Dy, r, boxes, box_limit, boxes_per_pixel_limit)

MAX_PIXELS_TO_IGNORE = 2;
BOX_MEMORY_FACTOR = 0.2;
BOXJITTER = true; % jitters box's center to lie closer to the object boundary
VIS = 999;

if ~exist('box_limit', 'var') || (box_limit < 0); box_limit = inf; end;
if ~exist('boxes_per_pixel_limit', 'var'); boxes_per_pixel_limit = 3; end;
[rows, cols, ~] = size(mask_in); imsize = [rows, cols];
dsk = strel('disk', 5);

mask = round(mask_in);
mask_vals = unique(mask_in);
n_mask_vals = length(mask_vals);

if isempty(boxes); n_boxes = 0; else; n_boxes = size(boxes, 1); end
diam = 2 * r + 1;
    
yy = repmat((-r:r)', 1, diam);
xx = repmat((-r:r) , diam, 1);
box_dist = sqrt(yy .* yy + xx .* xx);
inv_box_dist = 1./ max(1.0, box_dist);

% build edge image to cover
dx_l = Dx * mask(:);
dy_l = Dy * mask(:);
dx_l = reshape(dx_l, imsize);
dy_l = reshape(dy_l, imsize);
d_l  = dx_l + dy_l;
ad_l = abs(d_l);

look_img = abs(d_l) > 0;
look = look_img;

%-----------------------------------------------------------------------------
% go through existing boxes and choose to keep or toss
%-----------------------------------------------------------------------------
counts = zeros(imsize);

% sortboxes on confidence
if ~isempty(boxes);
  [~, boxes_ind] = sort(cat(1, boxes.conf), 1, 'descend');
  boxes = boxes(boxes_ind);
end

valid_boxes = true(n_boxes, 1);
for k = 1:n_boxes;
  bb = boxes(k);
  
  y = round(boxes(k).y);
  x = round(boxes(k).x);
  r = round(boxes(k).r);
  diam = 2 * r + 1;
  z = zeros(diam, diam);
  
  % modifies box's location based on the edge image
  if BOXJITTER;
    ymin = max(1,    y - r);
    ymax = min(rows, y + r);
    xmin = max(1,    x - r);
    xmax = min(cols, x + r);
    
    ys = ymin:ymax;
    xs = xmin:xmax;
    
    bedge = ad_l(ys,xs) > 0;
    [bedgey, bedgex] = find(bedge);
    if ~isempty(bedgey) && ~isempty(bedgex);
      bdist = (bedgey - boxes(k).y) .^ 2 + (bedgex - boxes(k).x) .^ 2;
      [~, eidx] = min(bdist);
      y = ymin + bedgey(eidx) - 1;
      x = xmin + bedgex(eidx) - 1;
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
  
  %---------------------------------------------------------------------------
  % changes to the box
  %---------------------------------------------------------------------------
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
  end

  %---------------------------------------------------------------------------
  % keep or toss box criteria
  %---------------------------------------------------------------------------
  % if box outside of image limits, drop
  if (y < 1) || (y > rows) || (x < 1) || (x > cols);
    valid_boxes(k) = false; continue;
  end

  % if fg or bg < %10 of the box area or if the learning mask is too small, 
  % then lower confidence
  msk = boxes(k).fg_mask > 0.5;
  msk_learn = msk .* boxes(k).fg_prob .* boxes(k).gmm_learn_colour_mask;
  nmsk = sum(msk(:));
  nmsk_learn = sum(msk_learn(:));
  TOOFEW = 0.1 * npx;
  TOOMUCH = 0.9 * npx;
  if (nmsk < TOOFEW) || (nmsk > TOOMUCH) || (nmsk_learn < TOOFEW) ...
    || (nmsk_learn > TOOMUCH);
    boxes(k).conf = boxes(k).conf - BOX_MEMORY_FACTOR;
  end
  
  % if lay_occr is in the background layer (i.e. occlusion constraint isn't
  % maintained), then drop
  if (boxes(k).lay_occr == 0) || (boxes(k).lay_occr <= boxes(k).lay_occd);
    valid_boxes(k) = false; continue;
  end
  
  % if confidence in colour model fall too low, drop
  if boxes(k).conf_colour < 0.25;
    valid_boxes(k) = false; continue;
  end
 
  % if general local classifier window confidence too low, drop
  if boxes(k).conf < 0.10;
    valid_boxes(k) = false; continue;
  end

  % if too many boxes in one spot, drop
  if any(vec(counts(ys, xs))) >= boxes_per_pixel_limit;
    valid_boxes(k) = false; continue;
  end

  %---------------------------------------------------------------------------
  % if kept, edit the look image and counts
  %---------------------------------------------------------------------------
  look(ys, xs) = 0;
  counts(ys, xs) = counts(ys, xs) + 1;
end

if ~isempty(valid_boxes);
  boxes = boxes(valid_boxes);
  n_boxes = sum(valid_boxes);
end

%-----------------------------------------------------------------------------
% add new boxes
%-----------------------------------------------------------------------------
p_all = find(look > 0);

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

k = 0;
nboxes = 0;
while sum(check_p) > MAX_PIXELS_TO_IGNORE;
  if size(p_out, 1) >= (box_limit - n_boxes); break; end

  if ~isempty(pts_to_add);
    % extract a point and add to output list
    p = pts_to_add(1);
    pts_to_add = pts_to_add(2:end);

    % drop boxes that go out of bounds
    [py, px] = ind2sub(imsize, p);
    ys = ((py - r):(py + r))';
    xs = ((px - r):(px + r))';
    if ((ys(1) < 1) || (ys(end) > rows) || (xs(1) < 1) || (xs(end) > cols));
      ys = ((max(1, py - r)):(min(rows, py + r)))';
      xs = ((max(1, px - r)):(min(cols, px + r)))';
      look(ys, xs) = 0;
      check_p = look(p_all);
      continue;
    end

    %-------------------------------------------------------------------------
    % get more box info
    %-------------------------------------------------------------------------
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
    
    assert(occluder > occluded, sprintf('%s: occluder/occluded', mfilename));
    
    %-------------------------------------------------------------------------
    % make box and add it
    %-------------------------------------------------------------------------
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
    
    new_boxes(nboxes, 1) = bb_new;
    
    % find new points to make boxes for
    p_new = edge_inds(look(edge_inds) > 0);
    if ~isempty(p_new);
      p_new = p_new(1);
      pts_to_add = [pts_to_add; p_new(:)];
    end

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
    out_img = draw_boxes(mask_in, p_to_draw, r);
    fig(55); clf; imagesc(out_img);
    drawnow; pause(0.1);
  end
end

if nboxes == 0; new_boxes = []; end
end

%-----------------------------------------------------------------------------
% helper functions
%-----------------------------------------------------------------------------
function edge_inds = put_box(imsize, box, sz)
[y, x] = deal(box(1), box(2));

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
