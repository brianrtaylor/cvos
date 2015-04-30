%---------------------------------------------------------------------------
% how this function is called
%
% [uvb_cbf, uvf_cbf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob] = cross_bilateral_filter_step_test_occ( ...
%   occb_prob, occf_prob, uvb, uvb_filt, uvf, uvf_filt, occb_rev, occf_rev, uvb_rev, uvf_rev, I0, I1, I2, m, VIS)
%
% designed to do cross bilateral filter in occ regions and spit out better occ regions
% 
%---------------------------------------------------------------------------
function [uvb_cbf, uvf_cbf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob] = cross_bilateral_filter_step_test_occ( ...
  occb_prob, occf_prob, uvb, uvb_filt, uvf, uvf_filt, occb_rev, occf_rev, uvb_rev, uvf_rev, I0, I1, I2, m, VIS)

if ~exist('VIS', 'var'); VIS = 100; end;
[rows,cols,~] = size(I1); imsize = [rows, cols];

% % moseg (sun)
% SIGMA_R_DIVIDE = 15.0;
% SIGMA_R_MIN = 0.5;
% SIGMA_R_MAX = 4.0;
% OCCMETHOD = 'LAB_RES';
% W = 25;

% segtrack (g4v)
OCCMETHOD = 'RGB_RES';
SIGMA_R_DIVIDE = 15.0;
SIGMA_R_MIN = 0.25;
SIGMA_R_MAX = 4.0;
SMALL_UV_MAG = 1.0;
W = 20;

W_imsize = ceil(min(rows, cols) / 20.0);
W = min(W, W_imsize);

ones_img = ones(imsize);

if ~exist('m', 'var'); m = true; end; % model: empty is ayvaci, else g4v 

dsk = strel('disk', 2);

if isempty(occf_rev); occf_rev = zeros(imsize); end;
if isempty(occb_rev); occb_rev = zeros(imsize); end;

% TODO: use backward occlusions (b_r2) to say that the pixels where the
% residual value in uvb > residual of uvf are places we shouldn't alter

%-------------------------------------------------------------------
% compute some variables (safe flow, occ{b,f}_prob, ruv{b,f})
%-------------------------------------------------------------------
uvb_safe = double(uvb);
uvb_safe(isnan(uvb) | isinf(uvb)) = 0;
uvb_rev_safe = uvb_rev;
uvb_rev_safe(isnan(uvb_rev) | isinf(uvb_rev)) = 0;
ruvb = max(vec(uvb_rev_safe)) - min(vec(uvb_rev_safe));

uvb_avg = abs(uvb_safe);
uvb_avg(uvb_avg < SMALL_UV_MAG) = [];
uvb_avg = mean(uvb_avg(:));

% occb_prob = occ_res_to_prob(occb);
if isempty(m); b_change = occb_prob; else; b_change = ones_img; end

sigma_r_b = clip(ruvb / SIGMA_R_DIVIDE, SIGMA_R_MIN, SIGMA_R_MAX);
fprintf('sigma_r_b: %0.2f\n', sigma_r_b);

uvf_safe = double(uvf);
uvf_safe(isnan(uvf) | isinf(uvf)) = 0;
uvf_rev_safe = uvf_rev;
uvf_rev_safe(isnan(uvf_rev) | isinf(uvf_rev)) = 0;
ruvf = max(vec(uvf_rev_safe)) - min(vec(uvf_rev_safe));

uvf_avg = abs(uvf_safe);
uvf_avg(uvf_avg < SMALL_UV_MAG) = [];
uvf_avg = mean(uvf_avg(:));

if isempty(m); f_change = occf_prob; else; f_change = ones_img; end

sigma_r_f = clip(ruvf / SIGMA_R_DIVIDE, SIGMA_R_MIN, SIGMA_R_MAX);
fprintf('sigma_r_f: %0.2f\n', sigma_r_f);

%-------------------------------------------------------------------
% run cross bilateral filter flow again with better occlusions
%-------------------------------------------------------------------
%-------------------------------------------------------------------
% create interesting flows
%-------------------------------------------------------------------
occf_prob(isnan(occf_prob) | isinf(occf_prob)) = 0.0;
occb_prob(isnan(occb_prob) | isinf(occb_prob)) = 0.0;

f_change = imclose(occf_prob, dsk);
b_change = imclose(occb_prob, dsk);
ob = repmat(b_change, [1, 1, 2]);
of = repmat(f_change, [1, 1, 2]);
% fig(1); imagesc([b_change, f_change]);

uvb_fun = uvb .* (1 - ob) + uvb_filt .* ob; % mostly uvb, little uvb_cbf
uvf_fun = uvf .* (1 - of) + uvf_filt .* of; % mostly uvf, little uvf_cbf
uvb_fun_2 = uvb_filt .* (1 - of) + uvb .* of; % mostly uvb_cbf, little uvb
uvf_fun_2 = uvf_filt .* (1 - ob) + uvf .* ob; % mostly uvf_cbf, little uvf

%-------------------------------------------------------------------
% get occlusions from "decent" flow
%-------------------------------------------------------------------
occb_fun = occ_calc(I1, I0, uvb_fun, OCCMETHOD);
occf_fun = occ_calc(I1, I2, uvf_fun, OCCMETHOD);
occb_fun_prob = occ_res_to_prob(max(0.0, occb_fun - occf_fun));
occf_fun_prob = occ_res_to_prob(max(0.0, occf_fun - occb_fun));

occb_filt = occ_calc(I1, I0, uvb_filt, OCCMETHOD);
occf_filt = occ_calc(I1, I2, uvf_filt, OCCMETHOD);
occb_filt_prob = occ_res_to_prob(max(0.0, occb_filt - occf_filt - occb_rev));
occf_filt_prob = occ_res_to_prob(max(0.0, occf_filt - occb_filt - occf_rev));

occb_fun_2 = occ_calc(I1, I0, uvb_fun_2, OCCMETHOD);
occf_fun_2 = occ_calc(I1, I2, uvf_fun_2, OCCMETHOD);
occb_fun_prob_2 = occ_res_to_prob(max(0.0, occb_fun_2 - occf_fun_2));
occf_fun_prob_2 = occ_res_to_prob(max(0.0, occf_fun_2 - occb_fun_2));

if VIS < 150;
  fig(244); imagesc([[occb_filt_prob; occb_fun_prob; occb_fun_prob_2], ...
    [occf_filt_prob; occf_fun_prob; occf_fun_prob_2]]);
end

%-------------------------------------------------------------------
% cross bilateral filter that flow
%-------------------------------------------------------------------
if (isnan(uvb_avg) || isnan(uvf_avg));
  [uvb_cbf, uvf_cbf] = deal(uvb_fun_2, uvf_fun_2);
else
  atic = tic();
  uvb_cbf = cbf_flow(uvb_fun_2, uvf, b_change, ...
    occb_filt_prob, occf_filt_prob, W, sigma_r_b, sigma_r_f, 2*uvb_avg);
  atic = toc(atic); fprintf('uvb cbf 2 (%0.3f s)\n', atic);

  atic = tic();
  uvf_cbf = cbf_flow(uvf_fun_2, uvb, f_change, ...
    occf_filt_prob, occb_filt_prob, W, sigma_r_f, sigma_r_b, 2*uvf_avg);
  atic = toc(atic); fprintf('uvf cbf 2 (%0.3f s)\n', atic);
end

if VIS < 150;
  fig(243); imagex([[uvb, uvf];[uvb_filt, uvf_filt]; ...
    [uvb_fun_2, uvf_fun_2];[uvb_cbf, uvf_cbf]]);
  title('uvb / filt / fun-2 / cbf | uvf / filt / fun-2 / cbf '); colorbar; drawnow; 
end

%-------------------------------------------------------------------
% cross bilateral filter those occlusions
%-------------------------------------------------------------------
occb_cbf = occ_calc(I1, I0, uvb_cbf, OCCMETHOD);
occf_cbf = occ_calc(I1, I2, uvf_cbf, OCCMETHOD);
  occb_cbf_prob = occ_res_to_prob(occb_cbf, true);
  occf_cbf_prob = occ_res_to_prob(occf_cbf, true);

occb_0 = max(0.0, occb_cbf - occf_cbf - occf_rev);
occf_0 = max(0.0, occf_cbf - occb_cbf - occb_rev);

occb_cbf_sole = cbf_occ(occb_0, uvb_cbf, uvf_safe, uvb_rev_safe, ...
  im2double(I1), max(uvb_avg / 2, 3.0), uvb_avg, sigma_r_b, 25, 0.25);
occf_cbf_sole = cbf_occ(occf_0, uvf_cbf, uvb_safe, uvf_rev_safe, ...
  im2double(I1), max(uvf_avg / 2, 3.0), uvf_avg, sigma_r_f, 25, 0.25);

% standard what we did before
occb_cbf_prob_2  = occ_res_to_prob(max(0.0, occb_cbf_sole - occf_cbf_sole));
occf_cbf_prob_2  = occ_res_to_prob(max(0.0, occf_cbf_sole - occb_cbf_sole));
occb_cbf_prob_2  = max(0.0, occb_cbf_prob_2 - occf_cbf_prob_2);
occf_cbf_prob_2  = max(0.0, occf_cbf_prob_2 - occb_cbf_prob_2);

%-------------------------------------------------------------------
% sanitize output
%-------------------------------------------------------------------
uvf_cbf(isnan(uvf) | isinf(uvf)) = inf;
uvb_cbf(isnan(uvb) | isinf(uvb)) = inf;

nooccb = sum(isinf(uvb_cbf), 3) > 0;
nooccf = sum(isinf(uvf_cbf), 3) > 0;

occb_cbf_prob = occb_cbf_prob_2;
occf_cbf_prob = occf_cbf_prob_2;
occb_cbf_prob(isnan(occb_cbf_prob) | isinf(occb_cbf_prob) | nooccb) = 0;
occf_cbf_prob(isnan(occf_cbf_prob) | isinf(occf_cbf_prob) | nooccf) = 0;

occb_cbf(isnan(occb_cbf) | isinf(occb_cbf) | nooccb) = 0;
occf_cbf(isnan(occf_cbf) | isinf(occf_cbf) | nooccf) = 0;
end
