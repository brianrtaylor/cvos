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
%----------------------------------------------------------------------------%
% Processing
% 
% things I tried
% * medfilt result
%
% TODO: things I should try
% * using another norm instead of L2, maybe L1?
%----------------------------------------------------------------------------%

switch DETECT_OCC_METHOD;
  case 'RGB_RES'
    r = compute_residual(I0, I1, uvf);
    occ = sqrt(sum(r .^ 2, 3));

  case 'LAB_RES'
    r = compute_residual(vl_xyz2lab(vl_rgb2xyz(I0)), ...
      vl_xyz2lab(vl_rgb2xyz(I1)), uvf);
    occ = sqrt(sum(r .^ 2, 3));

  case 'LAB_RES_PATCH'  
    r = compute_residual_btay(vl_xyz2lab(vl_rgb2xyz(I0)), ...
      vl_xyz2lab(vl_rgb2xyz(I1)), uvf, 'patch');
  
    occ = sqrt(sum(r .^ 2, 3));
  
  case 'RGB_RES_HIST'
    r = compute_residual_btay(I0, I1, uvf, 'hist');
    r2 = sum(abs(r), 3); % L1 - maybe better for histograms

    r2(isnan(r2)) = 0;
    occ = r2;
    % occ = sqrt(sum(r .^ 2, 3)); % chi^2?

  case 'LAB_RES_HIST'  
    r = compute_residual_btay(vl_xyz2lab(vl_rgb2xyz(I0)), ...
        vl_xyz2lab(vl_rgb2xyz(I1)), uvf, 'hist');

    r2 = sum(abs(r), 3); % L1 - maybe better for histograms

    r2(isnan(r2)) = 0;
    occ = r2;
    % occ = sqrt(sum(r .^ 2, 3)); % chi^2?
    
  case 'LAB_RES_MAGIC'
    r = compute_residual(vl_xyz2lab(vl_rgb2xyz(I0)), ...
      vl_xyz2lab(vl_rgb2xyz(I1)), uvf);
% %       r2 = zeros(M, N, D); r3 = r2;
% %       for k = 1:D;
% %         r2(:,:,k) = filter_occlusions_with_divuv(r(:,:,k), uvf, 0.04);
% %         r3(:,:,k) = medfilt2(r2(:,:,k), [2, 2]);
% %       end
% %       
% %       r4 = sqrt(sum(r3 .^ 2, 3));
    r2 = sum(r .^ 2, 3);
    r3 = filter_occlusions_with_divuv(r2, uvf, 0.05);
    r4 = medfilt2(r3, [2, 2]);
    occ = r4;

  case 'SPARSE_OCCLUSION'
    occ = -1;

  case 'FUN'
    % rgb
    rrgb = compute_residual(I0, I1, uvf);
    rrgb_ = sqrt(sum(rrgb .^ 2, 3));

    rlab = compute_residual(vl_xyz2lab(vl_rgb2xyz(I0)), ...
      vl_xyz2lab(vl_rgb2xyz(I1), uvf));
    rlab_ = sqrt(sum(rlab .^ 2, 3));

    occ =  (rrgb_ + rlab_) / 2;

  otherwise
    fprintf('occdetect: ERROR: unknown occ detection method: %s\n', DETECT_OCC_METHOD);
end

% postprocessing for NaNs (disocclusions)
occ(isnan(occ)) = 0;

% % filter occlusions with divuv
% Rb = size(occb, 3); occb2 = zeros(size(occb)); 
% for t = 1:Rb;
%   occb2(:,:,t) = filter_occlusions_with_divuv(occb(:,:,t), uvb(:,:,:,t), -10); 
% end
% 
% Rf = size(occf, 3); occf2 = zeros(size(occf));
% for t = 1:Rf;
%   occf2(:,:,t) = filter_occlusions_with_divuv(occf(:,:,t), uvf(:,:,:,t), -10); 
% end

% startDisplayImage([makeCentered(occb, occf, 3), makeCentered(occb2, occf2, 3)]);
% occb = occb2; % occf = occf2;
end
