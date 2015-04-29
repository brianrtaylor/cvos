% constraints_ = propagate_occlusion_pairs_forward(prev_constraints, uvb_rev);
%   where:
%       prev_constraints -- defined on the domain of I(t-1)
%       uvb_rev -- FLOW( I(t-1), I(t) )
%       constraints -- will be defined in I(t)
% the warping is done by the motion of occluder only, and the occluded point's   
% location is determined by keeping it the same distance (and orientation)
% away from occluder, as in the previous frame.
%
function [constraints_, w_uv, valid] = warp_constraints_forward(constraints, uvb_rev, SAFESPEEDSQUARED)
if ~exist('SAFESPEEDSQUARED', 'var'); SAFESPEEDSQUARED = 10.0; end;
imsize = [size(uvb_rev,1), size(uvb_rev,2)];
ux = uvb_rev(:,:,1);
uy = uvb_rev(:,:,2);
% consider only the motion of the occluder:
[X, Y] = meshgrid(1:imsize(2), 1:imsize(1));
xy_occd = [X(constraints(:,2)), Y(constraints(:,2))];
xy_occr = [X(constraints(:,1)), Y(constraints(:,1))];   
% vector pointing from occluded to occluder
v = xy_occr - xy_occd;

% confidence(s)
% so at squared uv_mag_occr, the falloff weight is about 0.3 (exp(-1.0))
uv_mag_occr_squared = ux(constraints(:,1)).^2 + uy(constraints(:,1)).^2;
uv_mag_occr = sqrt(uv_mag_occr_squared);
uv_mag_occr_avg = mean(uv_mag_occr(:));
% w_uv = exp( -uv_mag_occr / max(sqrt(SAFESPEEDSQUARED), uv_mag_occr_avg) );
w_uv = exp( -uv_mag_occr / max(sqrt(SAFESPEEDSQUARED), (uv_mag_occr_avg / 2.0)) );

% warp occluder:
x1_ = X(constraints(:,1)) + ux(constraints(:,1));
y1_ = Y(constraints(:,1)) + uy(constraints(:,1));
x2_ = x1_ - v(:,1);
y2_ = y1_ - v(:,2);    
% round everything
x1_ = round(x1_); x2_ = round(x2_);
y1_ = round(y1_); y2_ = round(y2_);
%
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
