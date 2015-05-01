function [w, w_div] = make_pixelwise_weights( idx, I, uvb_, uvf_, gpb, weights, opts)

%--------------------------------------------------------------------
% settings 
%--------------------------------------------------------------------
DIVFLAG = false;
FERRARIFLAG = true; 
LOCALFLOWADJUST = false; % not great (should try with more than cube too)
LOCALFLOWWEIGHTADJUST = false;
WEIGHT_MAG_THRESHOLD = 0.65;
THETA = 0.2;
MINWUV = 5e-3;
v2struct(opts);

weights = weights/sum(weights);

[rows, cols, chan] = size(I);
imsize = [rows, cols];
n = rows * cols;

%--------------------------------------------------------------------
% meat
%--------------------------------------------------------------------
I = reshape(I, [numel(I)/chan, chan]);
uvf = reshape(uvf_, [numel(uvf_)/2, 2]);
uvb = reshape(uvb_, [numel(uvb_)/2, 2]);
idx1 = idx(:,1);
idx2 = idx(:,2);

% rgb intensity weight
if weights(1)>0
  w_img = sum( (I(idx1,:) - I(idx2,:)).^2, 2);
  w_img = exp( -w_img/mean(w_img) );
else
  w_img = zeros( length(idx1), 1);
end

w_uvadjust = 0.0;
if weights(2)>0
  w_uvf = sum( (uvf(idx1,:) - uvf(idx2,:)).^2, 2);
  w_uvb = sum( (uvb(idx1,:) - uvb(idx2,:)).^2, 2);
  if ~LOCALFLOWADJUST;
    % flow weight (forward)
    mean_w_uvf = mean(w_uvf);
    fprintf('mean uvf: %0.3f\n', mean_w_uvf);
    w_uvf_ = exp( -w_uvf/ max(MINWUV, mean(w_uvf)) );

    % flow weight (backward)
    mean_w_uvb = mean(w_uvb);
    fprintf('mean uvb: %0.3f\n', mean_w_uvb);
    w_uvb_ = exp( -w_uvb/ max(MINWUV, mean(w_uvb)) );

    w_mag = 0.5*(w_uvf_+w_uvb_);
  else
    sz = [15, 15];
    mf = fspecial('average', sz);
    med_uvf = imfilter(reshape(w_uvf, [rows, cols]), mf);
    wmed_uvf = mean([med_uvf(idx1), med_uvf(idx2)], 2);
    w_uvf_ = exp( -w_uvf ./ max(MINWUV, wmed_uvf) );

    med_uvb = imfilter(reshape(w_uvb, [rows, cols]), mf);
    % med_uvb = medfilt2(reshape(w_uvb, [rows, cols]), sz);
    wmed_uvb = mean([med_uvb(idx1), med_uvb(idx2)], 2);
    w_uvb_ = exp( -w_uvb ./ max(MINWUV, wmed_uvb) );

    w_mag = 0.5*(w_uvf_ + w_uvb_);
  end

  % compute angle weights:
  if FERRARIFLAG;
    df = flow_angle_distance_mex( reshape(uvf, [rows, cols, 2]) );
    db = flow_angle_distance_mex( reshape(uvb, [rows, cols, 2]) );
    theta_dist = (df + db)/2;
    w_theta = exp( -theta_dist(:)/THETA );
    w_theta = repmat(w_theta,[2,1]);
    w_uv = w_mag .* ( w_theta.*(w_mag <= WEIGHT_MAG_THRESHOLD) + (w_mag > WEIGHT_MAG_THRESHOLD) );
  else
    w_uv = w_mag;
  end

  if LOCALFLOWWEIGHTADJUST;
    % based on flow differences
    wadj_uvf = exp(- w_uvf / WARP_SAFESPEEDSQUARED);
    wadj_uvb = exp(- w_uvb / WARP_SAFESPEEDSQUARED);
    w_uvadjust_dif = 0.5*wadj_uvf+0.5*wadj_uvb;

    % based on flow magnitude
    uvf_mag = sum(uvf_ .^ 2, 3);
    uvb_mag = sum(uvb_ .^ 2, 3);
    wadj_uvf = exp(- uvf_mag / WARP_SAFESPEEDSQUARED);
    wadj_uvb = exp(- uvb_mag / WARP_SAFESPEEDSQUARED);
    w_uvadjust_mag = 0.5*wadj_uvf(:)+0.5*wadj_uvb(:);
    if numel(idx1) ~= numel(w_uvadjust_mag)
      w_uvadjust_mag = repmat(w_uvadjust_mag,[2 1]);
    end
    
    % 0 = no change, 1 = remove all flow influence
    w_uvadjust = 0.5 * w_uvadjust_mag + 0.5 * w_uvadjust_dif;
    w_uvadjust = MAXUVADJUSTPROPORTION * w_uvadjust(:);
  end
else
  w_uv = zeros( length(idx1), 1);        
end
% keyboard;

% gpb weight:
% ???
if weights(3) > 0
  gpb = gpb(:);
  % w_gpb = abs(gpb(idx1)-gpb(idx2));
  % w_gpb = min(gpb(idx1),gpb(idx2));
  w_gpb = gpb(idx1);
else
  w_gpb = zeros( length(idx1), 1);
end
% w = w_img*weights(1) + (w_uvf+w_uvb)*weights(2)/2 + w_gpb*weights(3);
w_old = w_img*weights(1) + w_uv*weights(2) + w_gpb*weights(3);
w_wrong = 3*0.5*(1 - w_uvadjust) .* w_img*weights(1) ...
  + 3*w_uvadjust .* w_uv*weights(2) + 3*0.5*(1 - w_uvadjust) .* w_gpb*weights(3);

weights_adjusted = cat(2, ...
  (1 + w_uvadjust) * weights(1), ...
  (1 - w_uvadjust) .* weights(2), ...
  (1 + w_uvadjust) * weights(3));

weights_adjusted = weights_adjusted ./ repmat(sum(weights_adjusted, 2), 1, 3);

w = w_img .* (weights_adjusted(:,1)) ...
  + w_uv  .* (weights_adjusted(:,2)) ...
  + w_gpb .* (weights_adjusted(:,3));

%-------------------------------------------------------------------
% divergence scaling factor
% 
% TODO: I shouldn't be taking this max, but really associating
% w_divf with constraints from the forward flow, and w_divb with
% constraints from backward flow, and using those to decide. else
% w_divb can be high for an incorrect constraint from uvf (+vice versa)
%-------------------------------------------------------------------
w_div = zeros([rows, cols, 2]);
if DIVFLAG;
  [y1, x1] = ind2sub([rows, cols], idx1);
  [y2, x2] = ind2sub([rows, cols], idx2);
  
  vector_flow_uvf = uvf(idx1,:) - uvf(idx2,:);
  vector_flow_uvb = uvb(idx1,:) - uvb(idx2,:);
  vector_loc = ([x1, y1] - [x2, y2]);
  vector_loc_denom = repmat(sqrt(sum(vector_loc .^ 2, 2)), [1, 2]);
  vector_loc = vector_loc ./ vector_loc_denom; % L2 normalization
  
  divf = sum(vector_flow_uvf .* vector_loc, 2);
  divb = sum(vector_flow_uvb .* vector_loc, 2);

  DIVDENOM = 0.5;
  w_divf = max(0, 1 - exp( divf / DIVDENOM ));
  w_divb = max(0, 1 - exp( divb / DIVDENOM ));

  w_div = w_divf;
  % w_div = max(w_divf, w_divb);
  
  % TODO: should I incorporate distance to fall off with locations far away?
end

if VIS < 150;
  % debugging visuals if need be
  vis_w_img = reshape(min(w_img(1:n), w_img(n+1:2*n)), imsize);
  vis_w_uv  = reshape(min(w_uv(1:n) , w_uv(n+1:2*n)) , imsize);
  vis_w_gpb = reshape(min(w_gpb(1:n), w_gpb(n+1:2*n)), imsize);
  vis_w_old = reshape(min(w_old(1:n), w_old(n+1:2*n)), imsize);
  vis_w     = reshape(min(w(1:n)    , w(n+1:2*n))    , imsize);
  
  if LOCALFLOWWEIGHTADJUST;
    vis_w_uvadjust = reshape(w_uvadjust(1:n), imsize);
  end

  fig(2002); clf;
  vl_tightsubplot(2,3,1); imagesc(vis_w_img); notick; title('lab');
  vl_tightsubplot(2,3,2); imagesc(vis_w_uv); notick; title('uv');
  vl_tightsubplot(2,3,3); imagesc(vis_w_gpb); notick; title('gpb');

  vl_tightsubplot(2,3,4); imagesc(vis_w_old); notick;
  vl_tightsubplot(2,3,6); imagesc(vis_w); notick;
  vl_tightsubplot(2,3,5);
  if LOCALFLOWWEIGHTADJUST;
    imagesc(vis_w_uvadjust); notick;
  end
end

% OUTPUT
w = reshape(w, [rows, cols, 2]);
w_div = reshape(w_div, [rows, cols, 2]);
end
