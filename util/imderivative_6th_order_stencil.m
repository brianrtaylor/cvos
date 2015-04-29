function Ie = imderivative_6th_order_stencil(image, direction)

h = [-1, 9, -45, 0, 45, -9, 1] / 60;
if direction == 'y'
  h = h';
end

Ie = imfilter(image, h, 'corr', 'symmetric', 'same');
