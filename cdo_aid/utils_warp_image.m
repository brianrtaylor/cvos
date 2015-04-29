function [Iwarped] = utils_warp_image(I, w)
	[M, N, D] = size(I);
	[x,y] = meshgrid(1:N,1:M); 
	I = double(I);
  w = double(w);
	Iwarped = zeros(size(I));
	for k = 1:D
  		Iwarped(:,:,k) = vl_imwbackward(I(:,:,k), x+w(:,:,1), y+w(:,:,2) );
    end
end