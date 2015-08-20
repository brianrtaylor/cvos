%-----------------------------------------------------------------------------
% compute_residual  
%
% computes image residual between two images given a warping from I0 to I1
%
% @return: residual: residual image
% @return: Iwarped: warped I1 back to I0 reference frame
% @param: I0, I1: input images
% @param: w: warping from I0 to I1
%-----------------------------------------------------------------------------
function [residual, Iwarped] = compute_residual(I0, I1, w)

[M, N, D] = size(I0);
[x,y] = meshgrid(1:N,1:M); 

I0 = double(I0);
I1 = double(I1);

Iwarped = zeros(M, N, D);
for k = 1:D
	Iwarped(:,:,k) = interp2(I1(:,:,k), x+w(:,:,1), y+w(:,:,2));
end

residual = I0 - Iwarped;
end
