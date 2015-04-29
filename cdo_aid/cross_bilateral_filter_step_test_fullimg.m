%---------------------------------------------------------------------------
% how this function is called
%
% [uvb_cbf, uvf_cbf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob] = ...
%   cross_bilateral_filter_step(occb, occf, uvb, uvf, uvb_rev, uvf_rev, I0, I1, I2, m);
% 
%---------------------------------------------------------------------------
function [uvb_cbf, uvf_cbf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob] = cross_bilateral_filter_step_test_fullimg( ...
  occb, occf, uvb, uvf, uvb_rev, uvf_rev, I0, I1, I2, m, VIS)

if ~exist('VIS', 'var'); VIS = 100; end;

[rows,cols,~] = size(I1);

% segtrack (g4v)
OCCMETHOD = 'RGB_RES';
SIGMA_R_DIVIDE = 15.0;
SIGMA_R_MIN = 0.25;
SIGMA_R_MAX = 4.0;
SMALL_UV_MAG = 1.0;
W = 20;

W_imsize = ceil(min(rows, cols) / 20.0);
W = min(W, W_imsize);

i1 = im2double(I1);

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

occb_prob = occ_res_to_prob(occb, false, true);
% if isempty(m); b_change = occb_prob; else; b_change = ones_img; end

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

occf_prob = occ_res_to_prob(occf, false, true);
% if isempty(m); f_change = occf_prob; else; f_change = ones_img; end

sigma_r_f = clip(ruvf / SIGMA_R_DIVIDE, SIGMA_R_MIN, SIGMA_R_MAX);
fprintf('sigma_r_f: %0.2f\n', sigma_r_f);

%-------------------------------------------------------------------
% run cross bilateral filter flow
%-------------------------------------------------------------------
if (isnan(uvb_avg) || isnan(uvf_avg));
  [uvb_cbf, uvf_cbf] = deal(uvb_safe, uvf_safe);
else
  EXP_MAX_ = 5.0;
  EXP_NUM_ = 2000;
  expi = 1:EXP_NUM_;
  exp_fxn = exp( - expi ./ (expi * EXP_MAX_));
  
  atic = tic();
  % b_change = ones_img;
  b_change = occb_prob;
  % uvb_cbf = cbf_flow(uvb_safe, uvf_safe, b_change, ...
  %   occb_prob, occf_prob, W, sigma_r_b, sigma_r_f, 2*uvb_avg);
  uvb_cbf = cbf_flow(uvb_safe, uvf_safe, b_change, ...
    occb_prob, occf_prob, W, sigma_r_b, sigma_r_f, 2*uvb_avg, exp_fxn);
  atic = toc(atic); fprintf('uvf cbf (%0.3f s)\n', atic);
  
  atic= tic();
  % f_change = ones_img;
  f_change = occf_prob;
  uvf_cbf = cbf_flow(uvf_safe, uvb_safe, f_change, ...
    occf_prob, occb_prob, W, sigma_r_f, sigma_r_b, 2*uvf_avg, exp_fxn);
  atic = toc(atic); fprintf('uvb cbf (%0.3f s)\n', atic);
end

% fig(3); imagex([[uvb; uvb_cbf_ay; uvb_cbf], [uvf; uvf_cbf_ay; uvf_cbf]]);
% keyboard;

uvb_cbf_safe = uvb_cbf; uvb_cbf_safe(isinf(uvb_cbf) | isnan(uvb_cbf)) = 0;
uvf_cbf_safe = uvf_cbf; uvf_cbf_safe(isinf(uvf_cbf) | isnan(uvf_cbf)) = 0;

%-------------------------------------------------------------------
% run cross bilateral filter occlusions
%-------------------------------------------------------------------
occb_cbf = occ_calc(I1, I0, uvb_cbf, OCCMETHOD);
occf_cbf = occ_calc(I1, I2, uvf_cbf, OCCMETHOD);
%   occb_cbf_prob = occ_res_to_prob(occb_cbf, true);
%   occf_cbf_prob = occ_res_to_prob(occf_cbf, true);

occb_cbf_sole = cbf_occ(max(0.0, occb_cbf - occf_cbf), ...
  uvb_cbf_safe, uvf_safe, uvb_rev_safe, i1, ...
  max(uvb_avg * 2.0, 3.0), uvb_avg * 2.0, sigma_r_b, 25, 0.25);
% note: uvf_safe instead of uvf_cbf_safe - is this good?
occf_cbf_sole = cbf_occ(max(0.0, occf_cbf - occb_cbf), ...
  uvf_cbf_safe, uvb_safe, uvf_rev_safe, i1, ...
  max(uvf_avg * 2.0, 3.0), uvf_avg * 2.0, sigma_r_f, 25, 0.25);
% note: uvb_safe instead of uvb_cbf_safe - is this good?

occb_cbf_prob_v0 = occ_res_to_prob(occb_cbf_sole);
occf_cbf_prob_v0 = occ_res_to_prob(occf_cbf_sole);

%   occb_cbf_prob_ = occb_cbf_prob_v0 .* ( 1 - occf_cbf_prob_v0 );
%   occf_cbf_prob_ = occf_cbf_prob_v0 .* ( 1 - occb_cbf_prob_v0 );

%-------------------------------------------------------------------
% visualizations
%-------------------------------------------------------------------
if VIS < 250;
  fig(241); imagex([[uvb; uvb_cbf], [uvf; uvf_cbf]]);
  title('uvb / uvb-cbf | uvf / uvf-cbf');
  %     fig(242); imagesc([[occb_cbf_prob; occb_cbf_prob_v0; occb_cbf_prob_], ...
  %       [occf_cbf_prob; occf_cbf_prob_v0; occf_cbf_prob_]]);
  %     title(['occb-cbf-prob / occb-cbf-prob-v0 / occb-cbf-prob-' ...
  %       '| occf-cbf-prob / occf-cbf-prob-v0 / occf-cbf-prob-']);
  colorbar; drawnow;
end

%-------------------------------------------------------------------
% sanitize output
%-------------------------------------------------------------------
occb_cbf_prob = occb_cbf_prob_v0;
occf_cbf_prob = occf_cbf_prob_v0;
occb_cbf_prob(isnan(occb_cbf_prob) | isinf(occb_cbf_prob)) = 0;
occf_cbf_prob(isnan(occf_cbf_prob) | isinf(occf_cbf_prob)) = 0;
end
