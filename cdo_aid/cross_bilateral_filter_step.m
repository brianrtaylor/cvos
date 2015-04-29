%---------------------------------------------------------------------------
% how it's called
%
% [uvb_cbf, uvf_cbf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob] = ...
%   cross_bilateral_filter_step(occb, occf, uvb, uvf, uvb_rev, uvf_rev, I0, I1, I2, m);
% 
%---------------------------------------------------------------------------
function [uvb_cbf, uvf_cbf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob] = cross_bilateral_filter_step( ...
  occb, occf, uvb, uvf, uvb_rev, uvf_rev, I0, I1, I2, m)

  % moseg (ayvaci)
  SIGMA_R_DIVIDE = 15.0;
  SIGMA_R_MIN = 0.5;
  SIGMA_R_MAX = 4.0;
  W = 25;

  ones_img = ones(size(occf));

  if ~exist('m', 'var'); m = true; end; % model: empty is ayvaci, else g4v 

  % TODO: use backward occlusions (b_r2) to say that the pixels where the
  % residual value in uvb > residual of uvf are places we shouldn't alter

  %-------------------------------------------------------------------
  % compute some variables
  %-------------------------------------------------------------------
  uvf_safe = double(uvf);
  uvf_safe(isnan(uvf) | isinf(uvf)) = 0;
  uvf_rev_safe = uvf_rev;
  uvf_rev_safe(isnan(uvf_rev) | isinf(uvf_rev)) = 0;
  ruvf = max(vec(uvf_rev_safe)) - min(vec(uvf_rev_safe));

  occf_prob = occ_res_to_prob(occf);
  if isempty(m);
    f_change = occf_prob;
  else
    f_change = ones_img;
  end

  sigma_r_f = clip(ruvf / SIGMA_R_DIVIDE, SIGMA_R_MIN, SIGMA_R_MAX);
  fprintf('sigma_r_f: %0.2f\n', sigma_r_f);

  uvb_safe = double(uvb);
  uvb_safe(isnan(uvb) | isinf(uvb)) = 0;
  uvb_rev_safe = uvb_rev;
  uvb_rev_safe(isnan(uvb_rev) | isinf(uvb_rev)) = 0;
  ruvb = max(vec(uvb_rev_safe)) - min(vec(uvb_rev_safe));

  occb_prob = occ_res_to_prob(occb);
  if isempty(m);
    b_change = occb_prob;
  else
    b_change = ones_img; 
  end

  sigma_r_b = clip(ruvb / SIGMA_R_DIVIDE, SIGMA_R_MIN, SIGMA_R_MAX);
  fprintf('sigma_r_b: %0.2f\n', sigma_r_b);

  %-------------------------------------------------------------------
  % run cross bilateral filter flow
  %-------------------------------------------------------------------
  % spread / smooth occ
  uvb_avg = abs(uvb_safe);
  uvb_avg(uvb_avg < 1.0) = [];
  uvb_avg = mean(uvb_avg(:));
  % if isnan(uvb_avg); uvb_avg = 1.0; end;
  uvf_avg = abs(uvf_safe);
  uvf_avg(uvf_avg < 1.0) = [];
  uvf_avg = mean(uvf_avg(:));
  % if isnan(uvf_avg); uvf_avg = 1.0; end;

  if (isnan(uvb_avg) || isnan(uvf_avg));
    [uvb_cbf, uvf_cbf] = deal(uvb_safe, uvf_safe);
  else
    atic = tic();
    if isempty(m);
      uvb_cbf = cbf_flow_ayvaci(uvb_safe, uvf_safe, b_change, ...
        occb_prob, occf_prob, W, sigma_r_b, sigma_r_f, 2*uvb_avg);
    else
      uvb_cbf = cbf_flow(uvb_safe, uvf_safe, b_change, ...
        occb_prob, occf_prob, W, sigma_r_b, sigma_r_f, 2*uvb_avg);
    end
    atic = toc(atic); fprintf('uvf cbf (%0.3f s)\n', atic);
  
    atic = tic();
    if isempty(m);
      [uvf_cbf, g] = cbf_flow_ayvaci(uvf_safe, uvb_safe, f_change, ...
        occf_prob, occb_prob, W, sigma_r_f, sigma_r_b, 2*uvf_avg);
    else
      [uvf_cbf, g] = cbf_flow(uvf_safe, uvb_safe, f_change, ...
        occf_prob, occb_prob, W, sigma_r_f, sigma_r_b, 2*uvf_avg);
    end
    atic = toc(atic); fprintf('uvb cbf (%0.3f s)\n', atic);
  end

  %-------------------------------------------------------------------
  % run cross bilateral filter occlusions
  %-------------------------------------------------------------------
  occb_cbf = occ_calc(I1, I0, uvb_cbf, 'LAB_RES');
  uvb_cbf_safe = uvb_cbf;
  uvb_cbf_safe(isinf(uvb_cbf) | isnan(uvb_cbf)) = 0;

  occf_cbf = occ_calc(I1, I2, uvf_cbf, 'LAB_RES');
  uvf_cbf_safe = uvf_cbf;
  uvf_cbf_safe(isinf(uvf_cbf) | isnan(uvf_cbf)) = 0;

  occb_cbf_prob_0 = occ_res_to_prob(occb_cbf, true);
  occf_cbf_prob_0 = occ_res_to_prob(occf_cbf, true);

  % spread / smooth occ
  uvb_avg = abs(uvb_cbf_safe);
  uvb_avg(uvb_avg < 1.0) = [];
  uvb_avg = mean(uvb_avg(:));
  uvf_avg = abs(uvf_cbf_safe);
  uvf_avg(uvf_avg < 1.0) = [];
  uvf_avg = mean(uvf_avg(:));

  % TODO: use occf_rev, and opposite occb to discount occlusions found from occf
  % TODO: use occb_rev, and opposite occf to discount occlusions found from occb
  occb_cbf_sole = cbf_occ(occb_cbf, uvb_cbf_safe, uvf_safe, uvb_rev_safe, ...
    im2double(I1), max(uvb_avg / 2, 3.0), uvb_avg, sigma_r_b, max(50, mean(occf_cbf(:))), 0.25);
  % note: uvf_safe instead of uvf_cbf_safe - is this good?

  occf_cbf_sole = cbf_occ(occf_cbf, uvf_cbf_safe, uvb_safe, uvf_rev_safe, ...
    im2double(I1), max(uvf_avg / 2, 3.0), uvf_avg, sigma_r_f, max(50, mean(occb_cbf(:))), 0.25);
  % note: uvb_safe instead of uvb_cbf_safe - is this good?

  %-----------------
  % more fun / tries
  %%% occb_cbf_sole_v0 = cbf_occ(max(0.0, occb_cbf - occf_cbf), uvb_cbf_safe, uvf_safe, uvb_rev_safe, ...
  %%%   im2double(I1), max(uvb_avg / 2, 3.0), uvb_avg, sigma_r_b, max(50, mean(occf_cbf(:))), 0.25);
  %%% % note: uvf_safe instead of uvf_cbf_safe - is this good?

  %%% occf_cbf_sole_v0 = cbf_occ(max(0.0, occf_cbf - occb_cbf), uvf_cbf_safe, uvb_safe, uvf_rev_safe, ...
  %%%   im2double(I1), max(uvf_avg / 2, 3.0), uvf_avg, sigma_r_f, max(50, mean(occb_cbf(:))), 0.25);

  %%% occb_cbf_prob_v0 = occ_res_to_prob(occb_cbf_sole_v0);
  %%% occf_cbf_prob_v0 = occ_res_to_prob(occf_cbf_sole_v0);

  occb_cbf_prob_v0 = occ_res_to_prob(max(occb_cbf - occf_cbf, 0.0));
  occf_cbf_prob_v0 = occ_res_to_prob(max(occf_cbf - occb_cbf, 0.0));
  %-----------------
  

  % change to prob
  occb_cbf_prob = occ_res_to_prob(occb_cbf_sole);
  occf_cbf_prob = occ_res_to_prob(occf_cbf_sole);

  % fun with occlusions (using opposite occlusions to kill off some stuff)
  % TODO: use occf_rev, and opposite occb to discount occlusions found from occf
  % TODO: use occb_rev, and opposite occf to discount occlusions found from occb
  occb_cbf_prob_ = occ_res_to_prob(max(occb_cbf_sole - occf_cbf_sole, 0.0));
  occf_cbf_prob_ = occ_res_to_prob(max(occf_cbf_sole - occb_cbf_sole, 0.0));
  % occb_cbf_prob_ = max(occb_cbf_prob - occf_cbf_prob, 0.0);
  % occf_cbf_prob_ = max(occf_cbf_prob - occb_cbf_prob, 0.0);
  % occb_cbf_prob_ = occb_cbf_prob .* (1.0 - occf_cbf_prob);
  % occf_cbf_prob_ = occf_cbf_prob .* (1.0 - occb_cbf_prob);

  fig(240); imagesc([occb_cbf, occb_cbf_sole, occf_cbf, occf_cbf_sole]); 
  title('occb__cbf | occb__cbf__sole | occf__cbf | occf__cbf__sole'); colorbar; drawnow;
  fig(241); imagesc([occb_cbf_prob_0, occb_cbf_prob, occf_cbf_prob_0, occf_cbf_prob]);
  title('occb__cbf (gauss) | occb__cbf__sole (cbf occ) | occf__cbf (gauss) | occf__cbf__sole (cbf occ)');
  fig(242); imagesc([[occb_cbf_prob; occb_cbf_prob_; occb_cbf_prob_v0], ...
    [occf_cbf_prob; occf_cbf_prob_; occf_cbf_prob_v0]]);
  title(['occb__cbf__prob / occb__cbf__prob__ / occb__cbf__prob__v0 ' ...
    '| occf__cbf__prob / occf__cbf__prob__ / occf__cbf__prob__v0']);
  colorbar; drawnow;

  % occb_cbf_prob = occb_cbf_prob_v0;
  % occf_cbf_prob = occf_cbf_prob_v0;
  occb_cbf_prob = occb_cbf_prob_;
  occf_cbf_prob = occf_cbf_prob_;

  %%% %-------------------------------------------------------------------
  %%% % running second time with better flow and occlusion estimates 
  %%% %-------------------------------------------------------------------
  %%% %
  %%% % for however number of times you want to iterate it do the below
  %%% %
  %%% %-------------------------------------------------------------------
  %%% % run cross bilateral filter flow again with better occlusions
  %%% %-------------------------------------------------------------------
  %%% f_change_2 = occ_expand_mex(occf_cbf_prob, uvf_cbf);
  %%% atic = tic();
  %%% uvf_cbf_2 = cbf_flow(uvf_cbf, uvb, f_change_2, ...
  %%%   occf_cbf_prob, occb_cbf_prob, W, sigma_r_f, sigma_r_b, 2*uvf_avg);
  %%% atic = toc(atic); fprintf('uvb cbf (%0.3f s)\n', atic);

  %%% b_change_2 = occ_expand_mex(occb_cbf_prob, uvb_cbf);
  %%% atic = tic();
  %%% uvb_cbf_2 = cbf_flow(uvb_cbf, uvf, b_change_2, ...
  %%%   occb_cbf_prob, occf_cbf_prob, W, sigma_r_b, sigma_r_f, 2*uvb_avg);
  %%% atic = toc(atic); fprintf('uvf cbf (%0.3f s)\n', atic);

  %%% %-------------------------------------------------------------------
  %%% % run cross bilateral filter occlusions
  %%% %-------------------------------------------------------------------
  %%% occb_cbf_2 = occ_calc(I1, I0, uvb_cbf_2, 'LAB_RES');
  %%% occf_cbf_2 = occ_calc(I1, I2, uvf_cbf_2, 'LAB_RES');
  %%% occb_cbf_prob_0_2 = occ_res_to_prob(occb_cbf_2, true);
  %%% occf_cbf_prob_0_2 = occ_res_to_prob(occf_cbf_2, true);

  %%% % TODO: use occf_rev, and opposite occb to discount occlusions found from occf
  %%% % TODO: use occb_rev, and opposite occf to discount occlusions found from occb
  %%% occb_cbf_sole_2 = cbf_occ(occb_cbf_2, uvb_cbf_2, uvf_cbf_2, uvb_rev_safe, ...
  %%%   im2double(I1), max(uvb_avg / 2, 3.0), uvb_avg, sigma_r_b, max(50, mean(occf_cbf_2(:))), 0.25);
  %%% occf_cbf_sole_2 = cbf_occ(occf_cbf_2, uvf_cbf_2, uvb_cbf_2, uvf_rev_safe, ...
  %%%   im2double(I1), max(uvf_avg / 2, 3.0), uvf_avg, sigma_r_f, max(50, mean(occb_cbf_2(:))), 0.25);

  %%% % change to prob
  %%% occb_cbf_prob_2 = occ_res_to_prob(occb_cbf_sole_2);
  %%% occf_cbf_prob_2 = occ_res_to_prob(occf_cbf_sole_2);

  %%% occb_cbf_prob_2_ = occ_res_to_prob(max(occb_cbf_sole_2 - occf_cbf_sole_2, 0.0));
  %%% occf_cbf_prob_2_ = occ_res_to_prob(max(occf_cbf_sole_2 - occb_cbf_sole_2, 0.0));

  %%% fig(240); imagesc([occb_cbf, occb_cbf_sole, occb_cbf_sole_2, ...
  %%%   occf_cbf, occf_cbf_sole, occf_cbf_sole_2]); 
  %%% title(['occb__cbf | occb__cbf__sole | occb__cbf__sole__2 | occf__cbf '...
  %%%   '| occf__cbf__sole | occf__cbf__sole__2']); colorbar; drawnow;
  %%% fig(241); imagesc([occb_cbf_prob_0, occb_cbf_prob, occb_cbf_prob_2, ...
  %%%   occf_cbf_prob_0, occf_cbf_prob, occf_cbf_prob_2]);
  %%% title('???');
  %%% % title('occb__cbf (gauss) | occb__cbf__sole (cbf occ) | occf__cbf (gauss) | occf__cbf__sole (cbf occ)');
  %%% fig(242); imagesc([[occb_cbf_prob; occb_cbf_prob_; occb_cbf_prob_2_], ...
  %%%   [occf_cbf_prob; occf_cbf_prob_; occf_cbf_prob_2_]]);
  %%% title(['occb__cbf__prob / occb__cbf__prob__ / occb__cbf__prob__2__ | ' ...
  %%%   'occf__cbf__prob / occf__cbf__prob__ / occf__cbf__prob__2__']);
  %%% colorbar; drawnow;

  %%% occb_cbf_prob_2 = occb_cbf_prob_2_;
  %%% occf_cbf_prob_2 = occf_cbf_prob_2_;

  %-------------------------------------------------------------------
  % visuals and output
  %-------------------------------------------------------------------
  % fig(210); imagex([[uvb; uvb_cbf; uvb_cbf_2], [uvf; uvf_cbf; uvf_cbf_2]]); colorbar; 
  fig(210); imagex([[uvb; uvb_cbf], [uvf; uvf_cbf]]); colorbar; 
  title('uvb / uvb__cbf | uvf / uvf__cbf'); drawnow;
  % % fig(230); imagesc([[occb_prob, b_change, b_change_2];[occf_prob, f_change, f_change_2]]); colorbar;
  % fig(230); imagesc([[occb_prob, b_change];[occf_prob, f_change]]); colorbar;
  % title('occb__prob / b__change | occf__prob / f__change'); drawnow;

  uvf_cbf(isnan(uvf) | isinf(uvf)) = inf;
  uvb_cbf(isnan(uvb) | isinf(uvb)) = inf;
end
