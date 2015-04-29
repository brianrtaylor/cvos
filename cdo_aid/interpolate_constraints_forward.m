% constraints_ = propagate_occlusion_pairs_forward(prev_constraints, uvb_rev);
%   where:
%       prev_constraints -- defined on the domain of I(t-1)
%       uvb_rev -- FLOW( I(t-1), I(t) )
%       constraints -- will be defined in I(t)
% the warping is done by the motion of occluder only, and the occluded point's   
% location is determined by keeping it the same distance (and orientation)
% away from occluder, as in the previous frame.
%
function [constraints_, w_uv, valid] = interpolate_constraints_forward( ...
  constraints, uvb)
imsize = [size(uvb,1), size(uvb,2)];

% create occluder and occluded image
occdf = zeros(imsize);
occrf = zeros(imsize);
occdf(constraints(:,2)) = 1;
occrf(constraints(:,1)) = 1;

% warping begins
ux = uvb(:,:,1);
uy = uvb(:,:,2);
% consider only the motion of the occluder:
[X, Y] = meshgrid(1:imsize(2), 1:imsize(1));
xy_occd = [X(constraints(:,2)), Y(constraints(:,2))];
xy_occr = [X(constraints(:,1)), Y(constraints(:,1))];   
% vector pointing from occluded to occluder
v = xy_occr - xy_occd;

% warp the image
occrf_w = interp2(occrf, X + ux, Y + uy, 'nearest');
% occdf_w = interp2(occdf, X + ux, Y + uy); % this is wrong though, need to warp by occluder flow only :/
inds_occrf_w = find(occrf_w);
[y1, x1] = ind2sub(imsize, inds_occrf_w);
keyboard;
y2 = y1 + uy(constraints(:,1)) + v(:,2);
x2 = x1 + ux(constraints(:,1)) + v(:,1);
y2 = round(y2);
x2 = round(x2);
y1 = round(y1);
x1 = round(x1);


% confidence(s)
% so at squared uv_mag_occr, the falloff weight is about 0.3 (exp(-1.0))
SAFESPEEDSQUARED = 10.0;
uv_mag_occr_squared = ux(constraints(:,1)).^2 + uy(constraints(:,1)).^2;
w_uv = exp( -uv_mag_occr_squared / max(SAFESPEEDSQUARED, mean(uv_mag_occr_squared)) );

% check things
invalid = (x1_ < 1) | (x2_ < 1) | (y1_ < 1) | (y2_ < 1) | ...
          (x1_ > imsize(2)) | (x2_ > imsize(2)) | ...
          (y1_ > imsize(1)) | (y2_ > imsize(1));
valid = ~invalid;
%
x1_ = x1_(~invalid);
y1_ = y1_(~invalid);
x2_ = x2_(~invalid);
y2_ = y2_(~invalid);
%
idx1 = sub2ind(imsize, y1_, x1_);
idx2 = sub2ind(imsize, y2_, x2_);
%
constraints_ = [idx1(:), idx2(:)];
w_uv(invalid) = [];
end
