% I = imread('tulip.ppm');
load('/plot/vasiliy/CVPR15/data/moseg/cars7/flowsun/cars7_005.mat');
I = imread('/plot/vasiliy/CVPR15/data/moseg/cars7/cars7_05.ppm.png');

tic; out = recursive_bf_mex(I, 0.01, 0.1, 1, 5); toc;
figure(1);
imagesc([I, out]); axt;

tic; uvb_out = recursive_bf_mex( uvb, 0.05, 0.004, 0, 10); toc;
figure(2);
flowshow([uvb, uvb_out]);
