%----------------------------------------------------------------------------%
% 	COMPUTE_RESIDUAL  
% 		[RESIDUAL] = COMPUTE_RESIDUAL(I0, I1, W)
%----------------------------------------------------------------------------%
function [residual, Iwarped] = compute_residual(I0, I1, w)

[M, N, D] = size(I0);
[x,y] = meshgrid(1:N,1:M); 

I0 = double(I0); I1 = double(I1);

Iwarped = zeros(size(I0));
for k = 1:D
	Iwarped(:,:,k) = interp2(I1(:,:,k), x+w(:,:,1), y+w(:,:,2));
end

residual = I0 - Iwarped;

end
