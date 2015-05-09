%--------------------------------------------------------------------
% make_cvos_constraint_weights
%
% @param: opts : contains useful options
%   * relative_weights_uno_constraints : weights between intensity and flow 
%--------------------------------------------------------------------
function [constraints, constraint_weights, constraints_img] = make_cvos_unary_constraint_weights( ...
  constraints_in, constraint_weights_in, uvf, uvf_rev, I1, I2, opts)

% settings
DO_CONTRAINT_WEIGHT_SIGMOID = false;
relative_constraints_uno_weights = [0.5, 0.5];
v2struct(opts);

weights = relative_constraints_uno_weights ...
  / sum(relative_constraints_uno_weights);

% input
[rows, cols] = size(I1); 
imsize = [rows, cols];

% intensity change between I1, I1
% I1_ = reshape(I1, [numel(I1)/3, 3]);
% I2_ = reshape(I2, [numel(I2)/3, 3]);
% uvf_ = reshape(uvf, [numel(uvf)/2, 2]);
% uvf_rev_ = reshape(uvf_rev, [numel(uvf_rev)/2, 2]);
% idx = sub2ind(constrains_in, [rows, cols]);
idx = constraints_in;

if weights(1)>0
  w_img = sum( (I2(idx) - I1(idx)).^2, 2);
  w_img = exp( -w_img / mean(w_img) );
else
  w_img = zeros( length(idx), 1);
end

% flow sum magnitude between uvf, uvf_rev. if they don't = 0, that's a 
% stronger signal for occlusion
if weights(2) > 0;
  w_uv = (uvf(idx) + uvf_rev(idx)) .^ 2;
  w_uv = exp( -w_uv / mean(w_uv) );
else
  w_uv = zeros(length(idx),1);
end

constraint_weights_mod = w_img * weights(1) + w_uv * weights(2);
constraint_weights_mod = 1 - constraint_weights_mod;
constraint_weights = constraint_weights_mod .* constraint_weights_in;
constraints = constraints_in;

if DO_CONTRAINT_WEIGHT_SIGMOID;
  constraint_weights = sigmoid( ...
    constraint_weights, SIGMOID_MEAN, SIGMOID_SCALE);
end

% visuals and output
constraints_img = zeros(imsize);
% constraints_img_in = zeros(imsize);
% constraints_img_mod = zeros(imsize);

[constraints, constraint_weights] = aggregate_list( ...
  constraints, constraint_weights, false);
constraints_img(constraints)      = constraint_weights;
% constraints_img_in(constraints)   = constraint_weights_in;
% constraints_img_mod(constraints)  = constraint_weights_mod;

% fig(901); imagesc([constraints_img_in, constraints_img_mod, constraints_img]);
end
