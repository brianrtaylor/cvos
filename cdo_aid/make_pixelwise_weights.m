function w = make_pixelwise_weights( idx, I, uvb_, uvf_, gpb, weights, opts)

%--------------------------------------------------------------------
% settings 
%--------------------------------------------------------------------
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

if weights(2)>0
  w_uvf = sum( (uvf(idx1,:) - uvf(idx2,:)).^2, 2);
  w_uvb = sum( (uvb(idx1,:) - uvb(idx2,:)).^2, 2);

  % flow weight (forward)
  mean_w_uvf = mean(w_uvf);
  fprintf('mean uvf: %0.3f\n', mean_w_uvf);
  w_uvf_ = exp( -w_uvf/ max(MINWUV, mean(w_uvf)) );
  
  % flow weight (backward)
  mean_w_uvb = mean(w_uvb);
  fprintf('mean uvb: %0.3f\n', mean_w_uvb);
  w_uvb_ = exp( -w_uvb/ max(MINWUV, mean(w_uvb)) );
  
  w_mag = 0.5*(w_uvf_+w_uvb_);

  % compute angle weights
  df = flow_angle_distance_mex( reshape(uvf, [rows, cols, 2]) );
  db = flow_angle_distance_mex( reshape(uvb, [rows, cols, 2]) );
  theta_dist = (df + db)/2;
  w_theta = exp( -theta_dist(:)/THETA );
  w_theta = repmat(w_theta,[2,1]);
  w_uv = w_mag .* ( w_theta.*(w_mag <= WEIGHT_MAG_THRESHOLD) + (w_mag > WEIGHT_MAG_THRESHOLD) );
else
  w_uv = zeros( length(idx1), 1);        
end

% gpb weight:
if weights(3) > 0
  gpb = gpb(:);
  w_gpb = gpb(idx1);
else
  w_gpb = zeros( length(idx1), 1);
end

weights = weights ./ repmat(sum(weights, 2), 1, 3);
w = w_img .* (weights(:,1)) ...
  + w_uv  .* (weights(:,2)) ...
  + w_gpb .* (weights(:,3));

if VIS < 150;
  % debugging visuals if need be
  vis_w_img = reshape(min(w_img(1:n), w_img(n+1:2*n)), imsize);
  vis_w_uv  = reshape(min(w_uv(1:n) , w_uv(n+1:2*n)) , imsize);
  vis_w_gpb = reshape(min(w_gpb(1:n), w_gpb(n+1:2*n)), imsize);
  vis_w     = reshape(min(w(1:n)    , w(n+1:2*n))    , imsize);

  fig(2002); clf;
  vl_tightsubplot(2,3,1); imagesc(vis_w_img); notick; title('lab');
  vl_tightsubplot(2,3,2); imagesc(vis_w_uv); notick; title('uv');
  vl_tightsubplot(2,3,3); imagesc(vis_w_gpb); notick; title('gpb');

  vl_tightsubplot(2,3,6); imagesc(vis_w); notick;
  vl_tightsubplot(2,3,5);
end

% OUTPUT
w = reshape(w, [rows, cols, 2]);
end
