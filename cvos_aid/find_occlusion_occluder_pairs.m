function [occluded, occluder, outOfBounds] = find_occlusion_occluder_pairs(occ, uvf, uvb, MINDIST)

if ~exist('MINDIST', 'var'); MINDIST = 2.0; end;
if ~exist('VIS', 'var'); VIS = 100; end;

[M N] = size(occ);
mask = logical(occ > 0);
[X Y] = meshgrid(1:N, 1:M);
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
              | (isnan(sum([occluder, occluded], 2))) ... % removes NaNs
              | (sum(abs(occluder - occluded) < MINDIST, 2) == 2); % drops useless constraints

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

if VIS < 10;
fig(1010); imagesc([occ, occluded_img, -occluder_img, occluded_img - occluder_img]);
title('occ(red) | occd(red) | occr(blue) | occd(red), occr(blue)'); colorbar; drawnow;
end
end

%     uvb_warped = warp_quickly(uvb, uvf);
%     
%     uv = reshape(uvf+uvb_warped, [numel(uvf)/2, 2] );
%     
%     occluded = find(occ(:));
%     occD = zeros( length(occluded), 2);
%     occR = zeros( length(occluded), 2);
%     
%     for ii=1:length(occluded)
%         [yy xx] = ind2sub( size(occ), occluded(ii) );
%         occD(ii,:) = [yy,xx];
%         occR(ii,:) = [occD(ii,1)+uv(occluded(ii),1), occD(ii,2)+uv(occluded(ii),2)];
%     end
%     occR = round(occR);
%     invalid = ...
%         occR(:,1) <= 0 | occR(:,1) > size(occ,1) | ...
%         occR(:,2) <= 0 | occR(:,2) > size(occ,2);
%     occR = occR(~invalid,:);
%     occD = occD(~invalid,:);    
% end
% 
% function Iwarped = warp_quickly(I0, w)
% 	[M, N, D] = size(I0);
% 	[x,y] = meshgrid(1:N,1:M);
% 	I0 = double(I0);
% 
% 	Iwarped = zeros(size(I0));
% 	for k = 1:D
%   		Iwarped(:,:,k) = vl_imwbackward(I0(:,:,k), x+w(:,:,1), y+w(:,:,2) );
%     end
% end
