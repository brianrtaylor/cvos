function [dx, dy] = compute_directional_derivatives(in, sigma)

	% magic numbers
	sigma = 1;
	GaussianDieOff = .0001;  

	% Design the filters - a gaussian and its derivative
	pw = 1:30; % possible widths
	ssq = sigma^2;
	width = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff,1,'last');
	if isempty(width)
	  width = 1;  % the user entered a really small sigma
	end
	t = (-width:width);
	gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter

	% Find the directional derivative of 2D Gaussian (along X-axis)
	% Since the result is symmetric along X, we can get the derivative along
	% Y-axis simply by transposing the result for X direction.
	[x,y]=meshgrid(-width:width,-width:width);
	dgau2D=-x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);

	%smooth the image out
	aSmooth=imfilter(double(in),gau,'conv','replicate');   % run the filter across rows
	aSmooth=imfilter(aSmooth,gau','conv','replicate'); % and then across columns
	
	%apply directional derivatives
	dx = imfilter(aSmooth, dgau2D, 'conv','replicate');
	dy = imfilter(aSmooth, dgau2D', 'conv','replicate');
end
