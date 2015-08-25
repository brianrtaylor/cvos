%-----------------------------------------------------------------------------
% find_occlusion_occluder_pairs
%
% using occlusions, a warping, and it's reverse warping, this function yields
% the occluders for the occlusions, which make up the pairwise constraints
% for layer segmentation
%
% @return: occluded: valid occluded pixel indices
% @return: occluder: valid occluder pixel indices
% @return: outOfBounds: indices of constraints we removed
% @param: occ: occluded pixels
% @param: uvf: warping from I0 to Ix
% @param: uvb: warping from Ix to I0
% @param: MINDIST: minimum distance separating the occluded and occluder
%   pixels for a constraint to be valid (e.g. if an occluder pixels lands 
%   on the occluded pixel, then no useful constraint can be formed)
%-----------------------------------------------------------------------------
function [occluded, occluder, outOfBounds] = find_occlusion_occluder_pairs( ...
  occ, uvf, uvb, MINDIST)

if ~exist('MINDIST', 'var'); MINDIST = 2.0; end;
if ~exist('VIS', 'var'); VIS = 100; end;

[M, N] = size(occ);
mask = logical(occ > 0);
[X, Y] = meshgrid(1:N, 1:M);
uvf(isnan(uvf)) = inf;

% % warping attempt 1
wu = interp2(uvb(:,:,1), X + uvf(:,:,1), Y + uvf(:,:,2));
wv = interp2(uvb(:,:,2), X + uvf(:,:,1), Y + uvf(:,:,2));
wu(isnan(wu)) = inf; wv(isnan(wv)) = inf;
u = uvf(:,:,1) + wu;
v = uvf(:,:,2) + wv;

X0 = X + u;
Y0 = Y + v;

% constraints
occluder = [X0(mask), Y0(mask)];
occluded = [X(mask), Y(mask)];

occluder = round(occluder);
occluded = round(occluded);

outOfBounds =  ((occluder(:,1) < 1) | (occluder(:,1) > N) ...
              | (occluder(:,2) < 1) | (occluder(:,2) > M) ...
              | (occluded(:,1) < 1) | (occluded(:,1) > N) ...
              | (occluded(:,2) < 1) | (occluded(:,2) > M)) ...
              | (isnan(sum([occluder, occluded], 2))) ...
              | (sum(abs(occluder - occluded) < MINDIST, 2) == 2);

occluder = occluder(~outOfBounds,:);
occluded = occluded(~outOfBounds,:);
occluder = [occluder(:,2), occluder(:,1)];
occluded = [occluded(:,2), occluded(:,1)];

% visuals for show
occluded_img = zeros(M, N);
occluder_img = zeros(M, N);
occluded_ind = sub2ind([M, N], occluded(:,1), occluded(:,2));
occluder_ind = sub2ind([M, N], occluder(:,1), occluder(:,2));
occluded_img(occluded_ind) = 1;
occluder_img(occluder_ind) = 1;
end
