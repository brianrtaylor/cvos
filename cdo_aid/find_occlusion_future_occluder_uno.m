%--------------------------------------------------------------------
% find_occlusion_occluder_pairs
%
% finds pixels occluded in current frame that will become occluders
% in future frames, (unary occlusion term for detachable objects)
%--------------------------------------------------------------------
function [occluder, outOfBounds, occluder_img] = find_occlusion_future_occluder_uno(occ, uvf)

if ~exist('VIS', 'var'); VIS = 100; end;

[M N] = size(occ);
mask = logical(occ > 0);
[X Y] = meshgrid(1:N, 1:M);
uvf(isnan(uvf)) = inf;

% warping attempt 1
X0 = X + uvf(:,:,1);
Y0 = Y + uvf(:,:,2);

% constraints
occluder = [X0(mask), Y0(mask)];
occluder = round(occluder);

outOfBounds =  ((occluder(:,1) < 1) | (occluder(:,1) > N) ...
              | (occluder(:,2) < 1) | (occluder(:,2) > M) ...
              | (isnan(sum(occluder, 2)))); % removes NaNs

occluder = occluder(~outOfBounds,:);
occluder = [occluder(:,2), occluder(:,1)];

% visuals for show
occluder_img = zeros(M, N, 'single');
occluder_ind = sub2ind([M, N], occluder(:,1), occluder(:,2));
occluder_img(occluder_ind) = 1;
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
