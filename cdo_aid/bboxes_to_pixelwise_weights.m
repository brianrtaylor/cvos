%
% @return: ws : probability of a pixel after averaging all boxes over it
%   (2*npixels x 1 array)
% @return: cs : count of times a box touched a pixel (2*npixels x 1 array)
%
function [ws, cs, wxs, wys, cxs, cys] = bboxes_to_pixelwise_weights( ...
  imsize, boxes, FG, stub_conf, stub_conf_inter, stub_prob)

%----------------------------------------
% parse input
%----------------------------------------
n_boxes = size(boxes, 1);
wxs = zeros(imsize);
wys = zeros(imsize);
cxs = zeros(imsize);
cys = zeros(imsize);
rows = imsize(1);
cols = imsize(2);

if exist('stub_conf', 'var') && ~isempty(stub_conf);
  bb_conf = stub_conf;
else
  bb_conf = cat(1, boxes.conf_colour);
  % bb_conf = ones(n_boxes, 1);
end

if exist('stub_conf_inter', 'var') && ~isempty(stub_conf_inter);
  bb_conf_inter = stub_conf_inter;
else
  bb_conf_inter = cat(1, boxes.conf);
  % bb_conf_inter = ones(n_boxes, 1);
end

if ~exist('FG', 'var') || isempty(FG);
  FG = 1;
end

if exist('stub_prob', 'var') && ~isempty(stub_prob);
  bb_prob = stub_prob;
elseif FG == 1;
  bb_prob = cat(3, boxes.fg_prob);
else
  bb_prob = cat(3, boxes.bg_prob);
end

%----------------------------------------
% do work
%----------------------------------------
for k = 1:n_boxes;
  bb = boxes(k);
  
  %--------------------
  % gather indices
  %--------------------
  y = round(bb.y);
  x = round(bb.x);
  r = round(bb.r);
  
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

  ys1 = ys(1:(bheight - 1));
  % ys2 = ys(2:bheight);
  xs1 = xs(1:(bwidth - 1));
  % xs2 = xs(2:bwidth);
  
  bys1 = bys(1:(bheight - 1));
  bys2 = bys(2:bheight);
  bxs1 = bxs(1:(bwidth - 1));
  bxs2 = bxs(2:bwidth);

  %--------------------
  % edit weights
  %--------------------  
  pdiffx = abs(bb_prob(bys, bxs1, k) - bb_prob(bys, bxs2, k));
  pdiffy = abs(bb_prob(bys1, bxs, k) - bb_prob(bys2, bxs, k));
  wx = 1 - pdiffx;
  wy = 1 - pdiffy;
  
  if any(vec(isnan(wx))) || any(vec(isnan(wy)));
    fprintf('%s: NaNs inside of a box probability\n', mfilename); 
  end

  % combine via averaging
  new_wxs = wxs(ys, xs1) + bb_conf(k) * bb_conf_inter(k) * wx .* bb.invd(bys, bxs1);
  new_wys = wys(ys1, xs) + bb_conf(k) * bb_conf_inter(k) * wy .* bb.invd(bys1, bxs);
  
  if(any(vec(isnan(new_wxs))));
    fprintf('%s: NaNs inside of new_wxs\n', mfilename); 
  elseif (any(vec(isnan(new_wys))));
    fprintf('%s: NaNs inside of new_wys\n', mfilename); 
  end
  
  wxs(ys, xs1) = new_wxs;
  wys(ys1, xs) = new_wys;
 
  cxs(ys, xs1) = cxs(ys, xs1) + bb_conf(k) .* bb.invd(bys, bxs1);
  cys(ys1, xs) = cys(ys1, xs) + bb_conf(k) .* bb.invd(bys1, bxs);
end

%----------------------------------------
% output
%----------------------------------------
% wxs = wxs ./ max(1.0, cxs);
% wys = wys ./ max(1.0, cys);
wxs = wxs ./ cxs;
wys = wys ./ cys;
wxs(cxs == 0) = 0.0;
wys(cys == 0) = 0.0;
ws = cat(3, wxs, wys);
cs = cat(3, cxs, cys);
end
