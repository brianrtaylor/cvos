function [occout, div_uv_]  = filter_occlusions_with_divuv(occin, uv, thresh)

% FILTER_OCCLUSIONS_WITH_DIVUV filters the residual where div(uv) > -thresh
%  OCCOUT  = FILTER_OCCLUSIONS_WITH_DIVUV(OCCIN, UV, THRESH)
  div_uv = imderivative_6th_order_stencil(uv(:,:,1),'x') + ...
           imderivative_6th_order_stencil(uv(:,:,2),'y');
  h = fspecial('gaussian', [15 15], 6);
  div_uv_ = conv2(div_uv, h, 'same');

  occout = occin .* (div_uv_ < -thresh);

