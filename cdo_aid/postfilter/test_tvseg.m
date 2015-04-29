images = load_images('/home/vasiliy/data/kim_yu_na/', 'png', 1, inf);
load('/home/vasiliy/data/kim_yu_na/GT/GT.mat');

img = im2double(images{50});
gt = double( GT(:,:,50) );
gt(1:100,1:100) = 2;
gt(1:100,101:200) = 3;
unary = double( cat(3, gt==0, gt==1, gt==2, gt==3) );
%%
lambda = 2;
beta = .1;
max_iterations = 100;
img_gray = mean(img,3);
labels = bcv_tvsegment_mex( img_gray, -unary, lambda, beta, max_iterations);
% labels = 1-labels;

overlay1 = 0.3*img + 0.7*sc(double(gt),'jet');
overlay2 = 0.3*img + 0.7*sc(double(labels),'jet');

figure(1);
imagesc([overlay1, overlay2]); axt;

% figure(1);
% imagesc([overlay2]); axt;