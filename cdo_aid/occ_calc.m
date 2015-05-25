%----------------------------------------------------------------------------%
% detect_occ for single frame
%
% @param: I0: first image (warp to)
% @param: I1: second image (warp this image back to I0 ref frame) using uvf
% @param: params: parameters
% 
%----------------------------------------------------------------------------%
function [occ] = occ_calc(I0, I1, uvf, DETECT_OCC_METHOD)
if ~exist('DETECT_OCC_METHOD', 'var'); DETECT_OCC_METHOD = 'LAB_RES'; end

switch DETECT_OCC_METHOD;
  case 'RGB_RES'
    r = compute_residual(I0, I1, uvf);
    occ = sqrt(sum(r .^ 2, 3));

  case 'LAB_RES'
    r = compute_residual(vl_xyz2lab(vl_rgb2xyz(I0)), ...
      vl_xyz2lab(vl_rgb2xyz(I1)), uvf);
    occ = sqrt(sum(r .^ 2, 3));

  otherwise
    fprintf('occdetect: ERROR: unknown occ detection method: %s\n', DETECT_OCC_METHOD);
end

% postprocessing for NaNs (disocclusions)
occ(isnan(occ)) = 0;
end
