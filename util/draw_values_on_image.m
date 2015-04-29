%----------------------------------------------------------------------------%
% function I = draw_values_on_image(im, vals, alpha)
%
% draws the values (M x N double) on top of a grayscale version of the
% original image
%
% @param: im
% @param: vals
% @param: alpha : optional (mixing parameter)
%----------------------------------------------------------------------------%
function I = draw_values_on_image(im, vals, alpha)

if nargin < 3
	alpha = 0.5;
end

[M N D] = size(im);
MN = M * N;

MIN = double(min(min(vals)));
MAX = double(max(max(vals)));

jetcolors = jet;
nbins = size(jetcolors,1);
binsize = (MAX - MIN) / (nbins-1);

cmap = zeros([M N 3], 'uint8');
for k = 1:nbins
	binstart = MIN + binsize*(k-1)-binsize/2; 
	binend   = binstart + binsize;
	
	inds = find(vals >= binstart & vals < binend);
	c = jetcolors(k, :);
	cmap(inds)        = uint8(255*c(1));       % paint the region with c
	cmap(inds + MN)   = uint8(255*c(2));
	cmap(inds + 2*MN) = uint8(255*c(3));
end

if D > 1
  I = repmat(rgb2gray(im), [1 1 3]);
else
	I = repmat(im, [1 1 3]);
end

I = (1 - alpha) * I + alpha * cmap;
end
