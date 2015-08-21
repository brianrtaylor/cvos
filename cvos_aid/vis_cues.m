%-----------------------------------------------------------------------------
% vis_cues
%
% visualizes the constraints on top of the input image
%
% @return: W: output image with cues 
% @param: img (MxNx3): image to draw things on top of
% @param: constraints: occluder/occluded pairs to draw on image
% @param: constraint_weights: associated weights for each of the constraints
% @param: thresh: threshold to draw the constraints or not
%-----------------------------------------------------------------------------
function W = vis_cues(img, constraints, constraint_weights, thresh)
imsize = [size(img, 1), size(img, 2)];

if exist('thresh', 'var');
  keep = constraint_weights > thresh;
  constraints = constraints(keep, :);
  constraint_weights = constraint_weights(keep);
end

if isinteger(img); img = im2double(img); end;
if size(img, 3) == 1; 
  grayimg = repmat(img, [1 1 3]); 
  img = grayimg;
else
  grayimg = repmat(rgb2gray(img), [1 1 3]);
end

W = 0.3 * img + 0.7 * (grayimg / 2 + 0.5);

nconstraints = size(constraints, 1);
if ~exist('constraint_weights', 'var'); 
  constraint_weights = ones(nconstraints, 1);
end

if ~isempty(constraints);
  occd = zeros(imsize);
  occd(constraints(:,2))=constraint_weights;
  occd = cat(3, occd, zeros(imsize), zeros(imsize) ); % red
  
  occr = zeros(imsize);
  occr(constraints(:,1))=constraint_weights;
  occr = cat(3, occr, occr, zeros(imsize)); % yellow

  W = 0.1*W + 0.9*(W.*repmat(occd(:,:,1)==0,[1 1 3]) + occd);
  W = 0.1*W + 0.9*(W.*repmat(occr(:,:,1)==0,[1 1 3]) + occr);
end
end
