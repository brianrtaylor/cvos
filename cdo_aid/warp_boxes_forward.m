
% @param: boxes: struct with following fields
%   * center(k) = [y, x] : location of center in image
%   * radius : determines box size
%   * t(k) : struct with following fields
%     - mask : mask of foreground in the box
function [boxes_out, valid] = warp_boxes_forward( ...
  boxes, radii, layers, lay_occd, lay_occr, uvf, k)
if ~exist('k', 'var'); k = 1; end;
[rows, cols, ~] = size(uvf);
n = size(boxes, 1);
valid = true(n, 1);

boxes_out = [];
for k = 1:n;

  y = boxes(k, 1);
  x = boxes(k, 2);
  rad = radii(k);
  
  ymin = max(1,    round(y - rad));
  ymax = min(rows, round(y + rad));
  xmin = max(1,    round(x - rad));
  xmax = min(cols, round(x + rad));
  
  ysz = ymax - ymin + 1;
  xsz = xmax - xmin + 1;

  ys = ymin:ymax;
  xs = xmin:xmax;

  lay_fg = abs(layers(ys, xs) - lay_occr(k)) ...
    < abs(layers(ys, xs) - lay_occd(k));
  npx_box = (2 * rad + 1) ^ 2;
  if sum(lay_fg(:)) <= 0.02 * npx_box;
    lay_fg = true(ysz, xsz);
  end
  
  u = uvf(ys,xs,1);
  v = uvf(ys,xs,2);

  ub = mean(vec(u(lay_fg)));
  vb = mean(vec(v(lay_fg)));

  y_new = y + vb;
  x_new = x + ub;

  if (y_new >= 1) && (y_new <= rows) && (x_new >= 1) && (x_new <= cols);
    boxes_out = [boxes_out ; [y_new, x_new]];
  else
    valid(k) = false;
  end
end


% % % for k = 1:n;
% % %   b = boxes(k);
% % %   u = b.t(k).mask .* uvf(:,:,1);
% % %   v = b.t(k).mask .* uvf(:,:,2);
% % %   ub = mean(u(:));
% % %   vb = mean(v(:));
% % % 
% % %   b.center(k, :) = b.center(min(1, k-1),:) + [ub, vb];
% % %    
% % %   % remove boxes that go off the edge of the image
% % %   if (b.center(k, 1) < 1) || (b.center(k, 1) > rows) ...
% % %     || (b.center(k, 2) < 1) || (b.center(k, 2) > cols);
% % %     valid(k) = false;
% % %     boxes(k) = [];
% % %   end
% % % end
end
