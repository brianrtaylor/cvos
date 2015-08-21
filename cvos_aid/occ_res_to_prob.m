%-----------------------------------------------------------------------------
% occ_res_to_prob
%
% turns residual into a probability, after some smoothing
%
% @return: r2 (MxN): probability of occlusion
% @param: res (MxN): residual image
% @param: SMOOTH: whether to smooth the residual with a gaussian filter
% @param: DILATE: whether to dilate the residual via morphological op
%-----------------------------------------------------------------------------
function r2 = occ_res_to_prob(res, SMOOTH, DILATE)
% parameters and inputs
ALPHA = 2;
BETA = 10;
PLAY = false;
VIS = 100;
if ~exist('SMOOTH', 'var'); SMOOTH = false; end;
if ~exist('DILATE', 'var'); DILATE = false; end;

res2 = res;
res2(isnan(res2)) = 0.0;
if SMOOTH;
  g = fspecial('gauss', 5, 1);
  res2 = imfilter(res2, g);
end
if DILATE;
  g = fspecial('gauss', 3, 1);
  g = 0.5 * g / g(floor(numel(g) + 1) / 2);
  res2 = imfilter(res2, g);
end



% some playing
if PLAY;
  res_a1 = res;
  res_a1(abs(res_a1) < BETA) = [];
  BETA = mean(res_a1);
else
  rr = abs(res(:));
  
  max_res = max(rr);
  res_big = res(rr >= BETA);
  res_small = res(rr < BETA);
  if isempty(res_big);
    BETA = 2.0 * max_res / BETA;
  else
    BETA = max(BETA, mean(res_big(:)));
  end
end

r2 = 1 - exp(-(abs(res2 ./ BETA) .^ ALPHA) );
r2(isnan(res)) = 0;

% visuals
if VIS < 10;
  fig(1001); clf;
  vl_tightsubplot(1,3,1); 
  imagesc(res); colorbar; title('residual');
  vl_tightsubplot(1,3,2); 
  imagesc(res2); colorbar; title('residual smoothed');
  vl_tightsubplot(1,3,3); 
  imagesc(r2); colorbar; title('P(occ)');
end
end
