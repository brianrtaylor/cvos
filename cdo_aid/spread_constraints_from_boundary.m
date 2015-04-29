%
%
% function spreads constraints from boundary, so less chance of getting
% messed up after warp (warping across boundary, etc)
%
% @note: kind of a weaker perturb constraints, and only performed after the
% segmentation
% @note: perturb constraints should still be used otherwise
%
function constraints_out = spread_constraints_from_boundary( ...
  layers, constraints, rad)
constraints_out = constraints;

layers = round(layers);
[rows, cols, ~] = size(layers);
imsize = [rows, cols];

layer_edges = edge(double(layers), 'sobel');
g = fspecial('gauss', 2*rad+1, rad);
g = g / g(rad+1, rad+1);
dm_zone = imfilter(double(layer_edges), g) >= 1;

[y_occr, x_occr] = ind2sub(imsize, constraints(:, 1));
[y_occd, x_occd] = ind2sub(imsize, constraints(:, 2));

n = size(constraints, 1);
for k = 1:n;
  
  yr = round(y_occr(k));
  xr = round(x_occr(k));
  yd = round(y_occd(k));
  xd = round(x_occd(k));
  
  dy = round(yr - yd);
  dx = round(xr - xd);
  
  raddy = dy * rad ./ sqrt(dy * dy);
  raddx = dx * rad ./ sqrt(dx * dx);

  if dm_zone(yr, xr) > 0;
    % only go 1 length beyond initial location
    ey = yr + raddy;
    ex = xr + raddx;

    if abs(dy) < abs(dx);
      [oy, ox] = line_pts(yr, xr, ey, ex);
    else
      [ox, oy] = line_pts(xr, yr, ex, ey);
    end
    noy = length(oy);      
    
    ii = 1;
    while ii < noy && dm_zone(oy(ii), ox(ii)) > 0;
      ii = ii + 1;
    end
    
    constraints_out(k, 1) = sub2ind(imsize, oy(ii), ox(ii));
  end
  
  if dm_zone(yd, xd) > 0;
    % only go 1 length beyond initial location
    ey = yd - raddy;
    ex = xd - raddx;
    
    if abs(dy) <= abs(dx);
      [oy, ox] = line_pts(yd, xd, ey, ex);
    else
      [ox, oy] = line_pts(xd, yd, ex, ey);
    end
    noy = length(oy);
    
    ii = 1;
    while ii < noy && dm_zone(oy(ii), ox(ii)) > 0;
      ii = ii + 1;
    end
    
    constraints_out(k, 2) = sub2ind(imsize, oy(ii), ox(ii));
  end

end
end


% returns the first point starting from si, sy that returns a 0
function [oy, ox] = line_pts(sy, sx, ey, ex)
  
dy = ceil(double(ey) - double(sy));
dx = ceil(double(ex) - double(sx));

oy = zeros(dx, 1);
ox = zeros(dx, 1);

iy = 0;
ix = 0;

k = 0;

del = dy / dx;
err = 0;
inc = sign(dx);
for k = 1:abs(dx);
  ix = ix + inc;

  err = err + del;
  if err >= 0.5;
    iy = iy + inc;
    err = err - 1.0;
  end
  
  oy(k) = iy;
  ox(k) = ix;
end

oy = oy + sy;
ox = ox + sx;
end


% % % % % returns the first point starting from si, sy that returns a 0
% % % % function [oy, ox] = find_first_point_on_line(mask, sy, sx, ey, ex)
% % % % oy = round(sy);
% % % % ox = round(sx);
% % % % del = (double(ey) - double(sy)) / (double(ex) - double(sx));
% % % % err = 0;
% % % % while mask(sy, sx) > 0;
% % % %   ox = ox + 1;
% % % %   err = err + del;
% % % %   if err >= 0.5;
% % % %     oy = oy + 1;
% % % %     err = err - 1.0;
% % % %   end
% % % %   
% % % %   % got too far
% % % %   if oy > ey || ox > ex;
% % % %     oy = sy; ox = sx; return;
% % % %   end
% % % % end
% % % % end