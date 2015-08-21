%-----------------------------------------------------------------------------
% occ_calc
%
% detect occlusions for a single frame
%
% @return: occ (MxN): occlusion image
% @param: I0 (MxNx3): first image (warp to)
% @param: I1 (MxNx3): second image (warped back to I0 ref frame) using uvf
% @param: uvf (MxNx2): forward optical flow from I0 to I1
% @param: DETECT_OCC_METHOD: method for occlusion detection
% * LAB_RES: residual in lab colorspace
% * RGB_RES: residual in standard rgb colorspace
%-----------------------------------------------------------------------------
function occ = occ_calc(I0, I1, uvf, DETECT_OCC_METHOD)
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
    fprintf('occ_calc: unknown occ calc method: %s\n', DETECT_OCC_METHOD);
end

% postprocessing for NaNs (disocclusions)
occ(isnan(occ)) = 0;
end
