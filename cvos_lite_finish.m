%----------------------------------------------------------------------------%
% cdov_lite_finish(params, ADD, k)
%
% causal detachable objects for a single frame (i.e. first and last frames)
%
% @note: so far just does the last frame lol
%----------------------------------------------------------------------------%
function cdov_lite_finish(params, ADD, k)
setup;

%------------------------------------------------------------------
% setup
%------------------------------------------------------------------
img_path        = ADD.img_path;
flow_path       = ADD.flow_path;
outpath         = ADD.outpath;
seq             = ADD.seq;
files           = ADD.files;
flow_files      = ADD.flow_files;
nameStr         = ADD.nameStr;
edge_model      = ADD.edge_model;
past            = ADD.past;
boxes           = ADD.boxes;
dx_inds         = ADD.dx_inds;
dy_inds         = ADD.dy_inds;
Dx              = ADD.Dx;
Dy              = ADD.Dy;
objects         = ADD.objects;
object_map      = ADD.object_map;
problem         = ADD.problem;
FINISH          = ADD.FINISH;

object_mean_uvf_map     = ADD.object_mean_uvf_map;
% object_mean_uvf_map_t0  = ADD.object_mean_uvf_map_t0;

opts_flow_occ               = ADD.opts_flow_occ;
opts_current_weights        = ADD.opts_current_weights;
opts_constraint_weights     = ADD.opts_constraint_weights;
opts_propagate              = ADD.opts_propagate;
opts_unity                  = ADD.opts_unity;
opts_propagate_constraints  = ADD.opts_propagate_constraints;
opts_fg                     = ADD.opts_fg;
opts_boxes                  = ADD.opts_boxes;

fprintf('=========================== time: %d =============\n', k);

%------------------------------------------------------------------
% read all files around current image + pre-processing
%------------------------------------------------------------------
% kk = k; m1 = -1;
% if k > FINISH; kk = 1; m1 = 1; end

I1 = imread(fullfile(img_path, files(k).name));
I0 = imread(fullfile(img_path, files(k - 1).name));
if (size(I1,3)==1), I1 = repmat(I1,[1 1 3]); end;
if (size(I0,3)==1), I0 = repmat(I0,[1 1 3]); end;
[rows, cols, ~] = size(I1);
imsize = [rows, cols];
uvsize = [rows, cols, 2];
rgbsize = [rows, cols, 3];
i1 = im2double(rgb2gray(I1));
i0 = im2double(rgb2gray(I0));

I2 = zeros(rgbsize);
i2 = zeros(imsize);

% pre-processing on current t
I1_bflt = recursive_bf_mex(I1, 0.01, 0.1, 1, 5);
I0_bflt = recursive_bf_mex(I0, 0.01, 0.1, 1, 5);
i1_bflt = im2double(I1_bflt);
Ilab = vl_xyz2lab(vl_rgb2xyz( uint8(i1_bflt*255) ));

E = edgesDetect(I1_bflt, edge_model.model);
EdgeDollar = 0.75*(1-E);

%------------------------------------------------------------------
% some structures setup
%------------------------------------------------------------------
% past.uvf      = zeros(uvsize);
% past.uvf_rev  = zeros(uvsize);
% past.occf     = zeros(imsize);
% past.occf_rev = zeros(imsize);

unary_constraints_map = zeros(imsize);

%------------------------------------------------------------------
% tao_flow_occ: flow + occlusions + pre-processing
%------------------------------------------------------------------
[uvb, uvb_cbf, uvf, uvf_cbf, uvb_rev, uvf_rev, ...
  occb, occb_cbf, occb_cbf_prob, occf, occf_cbf, occf_cbf_prob, ...
  ~, ~, occb_rev_prob, occf_rev_prob] = ...
  tao_flow_occ(flow_path, flow_files, seq, k, FINISH, ...
  i0, i1, i2, I0, I1, I2, past, opts_flow_occ);

% % % if (k > FINISH);
% % %   tmp = uvf; uvf = uvb; uvb = tmp; % swap uvb, uvf
% % %   tmp = uvf_cbf; uvf_cbf = uvb_cbf; uvb_cbf = tmp;
% % %   tmp = uvf_rev; uvf_rev = uvb_rev; uvb_rev = tmp;
% % %   tmp = occf_cbf; occf_cbf = occb_cbf; occb_cbf = tmp;
% % %   tmp = occf_cbf_prob; occf_cbf_prob = occb_cbf_prob; occb_cbf_prob = tmp;
% % %   tmp = occf_rev_prob; occf_rev_prob = occb_rev_prob; occb_rev_prob = tmp; 
% % % end

uvf_bflt = double(uvf); % 0
uvf_cbf_bflt = double(uvf); % 0
uvb_bflt = recursive_bf_mex(double(uvb), 0.05, 0.004, 0, 10);
uvb_cbf_bflt = recursive_bf_mex(double(uvb_cbf), 0.05, 0.004, 0, 10);
uvb_rev_bflt = recursive_bf_mex(double(uvb_rev), 0.05, 0.004, 0, 10);

%------------------------------------------------------------------
% past.t0 edits (aka, warping things from the future to this frame)
%------------------------------------------------------------------
occb_mask = occb_cbf_prob > params.OCCPROBLAYERTHRESH;
% put into the frame 1 reference frame
layers_t0 = round(utils_warp_image(past.layers, uvb_cbf_bflt)); % warp back
[past.t0.layers, layers_t0_msfm] = remove_occb_mask( ...
  layers_t0, occb_mask, 0, 0); % past = t+1 lol

[boxes, prob_box_fg, prob_box_bg, count_box_fg, count_box_bg, ...
  weights_box] = tao_get_prob_from_boxes( ...
  occb_mask, past, i1_bflt, boxes, opts_boxes);

%------------------------------------------------------------------
% constraints now
%------------------------------------------------------------------
[constraints_now_b, constraint_weights_now_b, indbad_now_b] ...
  = make_tao_constraints(occb_cbf_prob, uvb_cbf_bflt, uvb_rev_bflt, ...
  params.OCCPROB, params.MINCONSTRAINTDIST, params.VIS, 16);
% % [constraints_now_f, constraint_weights_now_f, indbad_now_f] ...
% %   = make_tao_constraints(occf_cbf_prob, uvf_cbf_bflt, uvf_rev_bflt, ...
% %   params.OCCPROB, params.MINCONSTRAINTDIST, params.VIS, 17);
% [constraints_now_f, constraint_weights_now_f, indbad_now_f] = deal([], [], []);

% % if params.UNARYOCC;
% %   [constraints_uno_now_f, constraint_weights_uno_now_f, indbad_uno_now_f, ...
% %     constraints_uno_img_now_f] = make_tao_unary_constraints( ...
% %     occf_cbf_prob, uvf_cbf, params.OCCPROB);
% % end

%------------------------------------------------------------------
% weights now
%------------------------------------------------------------------
[weights_now, weights_inds] = make_tao_weights_current_frame( ...
  Ilab, uvb_cbf_bflt, uvf_cbf_bflt, dx_inds, dy_inds, ...
  EdgeDollar, opts_current_weights);

% % if params.BOXHELP && ~isempty(boxes);
% %   weights_now = (1 - params.BOXTVWEIGHT) * weights_now ...
% %     + max(0.0, params.BOXTVWEIGHT * weights_boxes);
% % end
  
%------------------------------------------------------------------
% constraints weights now
%
% TODO: replace I1_bflt with Ilab
%------------------------------------------------------------------
[constraints_now_b, constraint_weights_now_b, ...
  constraint_weights_now_nodiv_b] = make_tao_constraint_weights( ...
  constraints_now_b, constraint_weights_now_b, -uvb_cbf_bflt, uvb_cbf_bflt, ...
  I1_bflt, weights_now, EdgeDollar, opts_constraint_weights, 18);
% % [constraints_now_f, constraint_weights_now_f, ...
% %   constraint_weights_now_nodiv_f] = make_tao_constraint_weights( ...
% %   constraints_now_f, constraint_weights_now_f, uvb_bflt, uvf_cbf_bflt, ...
% %   I1_bflt, weights_now, EdgeDollar, opts_constraint_weights, 19);
[constraints_now_f, constraint_weights_now_f, ...
  constraint_weights_now_nodiv_f] = deal([], [], []);

% % if params.UNARYOCC;
% %   [constraints_uno_now_f, constraint_weights_uno_now_f, ...
% %     constraints_uno_img_now_f] = make_tao_unary_constraint_weights( ...
% %     constraints_uno_now_f, constraint_weights_uno_now_f, ...
% %     uvf, uvf_rev, I1, I2, opts_constraint_weights);
% % end
  
%------------------------------------------------------------------
% weights from prior frame (future)
%------------------------------------------------------------------
if ~isempty(past.weights) && params.CAUSAL && params.WEIGHTHELP;
  [weights, wx_l, wy_l, past] = propagate_tao_weights( ...
    weights_now, uvb_cbf_bflt, past, occb_mask, ...
    Dx, Dy, k, opts_propagate);
else
  weights = weights_now;
end
  
%---------------------------------------------------------------------
% Block that introduces a layer unity prior
%---------------------------------------------------------------------
if ~isempty(past.unity) && params.UNITYHELP && params.CAUSAL ...
  && params.PROB_UNITY > 0;
  unity = tao_unity_prior(past, uvb_bflt, uvb_cbf_bflt, ...
    wx_l, wy_l, occb_mask, k, opts_unity);
else
  unity = zeros([rows, cols, 2]);
end

%------------------------------------------------------------------
% constraints from prior frame
%
% TODO: pass Ilab in as well for the current frame weights computation
%------------------------------------------------------------------
if ~exist('prob_box_fg', 'var'); prob_box_fg = []; end;
if ~exist('count_box_fg', 'var'); count_box_fg = []; end;

if ~isempty(past.constraints_b) && params.CAUSAL && params.CONSTRAINTHELP;
  [constraints_b, constraint_weights_b, constraints_causal_b, ...
    constraint_weights_old_b, valid_b] = propagate_tao_constraint_weights( ...
    I0_bflt, I1_bflt, i1_bflt, uvb_cbf_bflt, -uvb_cbf_bflt, past, ...
    past.constraints_b, past.constraint_weights_b, past.weights, ...
    weights, past.xi_b, EdgeDollar, constraints_now_b, ...
    constraint_weights_now_b, opts_propagate_constraints, ...
    prob_box_fg, count_box_fg, objects, 11);
else
  constraints_b = constraints_now_b;
  constraint_weights_b = constraint_weights_now_b;
  constraints_causal_b = [];
  constraint_weights_old_b = [];
end

if ~isempty(past.constraints_f) && params.CAUSAL && params.CONSTRAINTHELP;
  [constraints_f, constraint_weights_f, constraints_causal_f, ...
    constraint_weights_old_f, valid_f] = propagate_tao_constraint_weights( ...
    I0_bflt, I1_bflt, i1_bflt, -uvb_bflt, uvb_bflt, past, ...
    past.constraints_f, past.constraint_weights_f, past.weights, ...
    weights, past.xi_f, EdgeDollar, constraints_now_f, ...
    constraint_weights_now_f, opts_propagate_constraints, ...
    [], [], objects, 12);
else
  constraints_f = constraints_now_f;
  constraint_weights_f = constraint_weights_now_f;
  constraints_causal_f = [];
  constraint_weights_old_f = [];
end

% more pruning based on duplicates
[constraints, constraint_weights, ~] = aggregate_pairs_fast( ...
  double([constraints_b; constraints_f]), ...
  double([constraint_weights_b; constraint_weights_f]));

%---------------------------------------------------------------------
%  Block that introduces a ``occlusion-driven foreground prior''
%---------------------------------------------------------------------
% TODO: constraint_weights_uno modifiying
if params.CAUSAL && params.UNARYOCC;
  [unary_constraints_map, past.unary_constraints_map] = ...
    make_unary_constraints(past, ...
    constraints_uno_img_now_b, constraints_uno_img_now_f, opts_unary);
  
  if params.VIS < 200;
    fig(451); imagesc(unary_constraints_map); colorbar; title(k);
  end
end

%---------------------------------------------------------------------
% Block that introduces a ``layer-driven foreground prior''
%---------------------------------------------------------------------
if ~isempty(past.prob_fg) && params.CAUSAL && params.PROB_FG > 0;
  prob_fg = tao_fg_prior(past, uvb_cbf_bflt, occb_mask, opts_fg);
else
  prob_fg = zeros(imsize);
end

%------------------------------------------------------------------
% solve problem
%------------------------------------------------------------------
if isempty(constraints);
  problem.constraints = [[1,2];[2,3]];
  problem.nocc_constraints = 2;
  problem.lambda = [0.0; 0.0];
else
  problem.nocc_constraints = size(constraints,1);
  problem.constraints = constraints;
  problem.lambda = params.LAMBDA * double(constraint_weights);
end
  
%------------------------------------------------------------------
% problem.Wx (weights)
%------------------------------------------------------------------
if params.CAUSAL && params.BOXHELP && ~isempty(boxes);
  weights_b4 = weights;
  weights = min(weights, max(0.0, ...
    -weights_box * params.BOXTVWEIGHT + double(weights_b4)));
  
  if params.UNITYHELP && ~isempty(unity);
    unity = min(unity, max(0.0, ...
      +weights_box * params.BOXTVWEIGHT + double(unity)));
  end
  
  if params.VIS < 150;
    ima = [-weights_box(:, :, 1), -weights_box(:, :, 2), max(-weights_box, [], 3)];
    imb = [weights_b4(:, :, 1), weights_b4(:, :, 2), max(weights_b4, [], 3)];
    imc = [weights(:, :, 1), weights(:, :, 2), max(weights, [], 3)];
    % box weights, weightsb4, weights+boxweights
    fig(9); clf; sc([ima; imb; imc], 'jet', [0, 2]); drawnow;
  end
end

weights_cut = max(0.0, 1.0 - max(0.0, double(weights - unity)));
weights_unity = max(0.0, double(unity - weights));
Wx = params.PAIR * (weights_cut + weights_unity);
problem.Wx = double(Wx(:));
problem.edges = weights_inds;

probwx_map = min(Wx, [], 3);
probwx_img = sc(clip(probwx_map, 0, 5), 'jet');
probwx_img = 0.6 * probwx_img + 0.4 * repmat(i1, [1,1,3]);

% visualization
if params.VIS < 160;
  fig(6); clf; imagesc(probwx_img);
  title(sprintf('problem.Wx (t = %d): %0.6f', k, max(Wx(:))));
end
if params.VIS < 110; % 210
  fig(7); clf; imagesc(probwx_map);
  title(sprintf('problem.Wx (t = %d): %0.6f', k, max(Wx(:))));
end
if params.VIS < 130; % 230 % visualization of combination here
  cut_img   = max(0.0, weights - unity);
  u_img     = max(0.0, unity - weights);
  w_cut_img = max(0.0, 1.0 - cut_img);
  w_img     = w_cut_img + u_img;
  
  im1 = [unity(:, :, 1)    , unity(:, :, 2)    , max(unity, [], 3)];
  im2 = [weights(:, :, 1)  , weights(:, :, 2)  , max(weights, [], 3)];
  im3 = [u_img(:, :, 1)    , u_img(:, :, 2)    , max(u_img, [], 3)];
  im4 = [cut_img(:, :, 1)  , cut_img(:, :, 2)  , max(cut_img, [], 3)];
  im5 = [w_cut_img(:, :, 1), w_cut_img(:, :, 2), min(w_cut_img, [], 3)];
  im6 = [w_img(:, :, 1)    , w_img(:, :, 2)    , min(w_img, [], 3)];
  
  % unity, edge_weight, unity-edge, edge-unity, cuts, final combo
  im_all = [[im1,im2];[im3,im4];[im5,im6]];
  fig(8); clf; sc(im_all, [0, 2], 'jet'); drawnow;
end

%------------------------------------------------------------------
% problem.kappa (fg prior)
%------------------------------------------------------------------
if params.CAUSAL && ~isempty(past.layers) ...
    && ((params.BOXHELP > 0) || (params.PROB_FG > 0));
  problem.USE_TEMPORAL_PENALTY = 1;
  
  if params.CAUSAL && params.BOXHELP && ~isempty(boxes);
    % use this, since the other will hurt 2 objects moving towards overlap
    box_kappa = params.PROB_BOX_FG * prob_box_fg;
    % box_kappa = params.PROB_BOX_FG * (prob_box_fg - prob_box_bg);
  else
    box_kappa = 0;
  end
  
  if params.UNARYOCC;
    unary_kappa = params.UNARY_LAMBDA * past.unary_constraints_map(:);
  else
    unary_kappa = 0;
  end
  
  fg_kappa = params.PROB_FG * prob_fg(:);
  
  problem.kappa = double(max(0.0, fg_kappa(:) + unary_kappa(:) + box_kappa(:)));
  
  if params.VIS < 300;
    kappa_img = sc(reshape(clip(problem.kappa, 0, 1), imsize), 'jet');
    kappa_img = 0.6 * kappa_img + 0.4 * repmat(i1, [1,1,3]);
    if params.BOXHELP && ~isempty(boxes);
      box_xy = bbox_to_box_center(boxes);
      box_conf = cat(1, boxes.conf);
      [~, box_img] = draw_boxes(past.t0.layers, box_xy, ...
        params.BOX_RAD, box_conf);
      b = find(box_img > 0);
      npx = rows*cols;
      kappa_img(b) = 0.5 * kappa_img(b) + 0.5;
      kappa_img(b + npx) = 0.5 * kappa_img(b + npx) + 0.5;
      kappa_img(b + 2*npx) = 0.5 * kappa_img(b + 2*npx) + 0.5;
    end
    fig(5); clf; imagesc(kappa_img);
    title(sprintf('problem.kappa (t = %d): %0.6f', ...
      k, max(vec(problem.kappa))));
    drawnow;
  end
else
  problem.USE_TEMPORAL_PENALTY = 0;
  problem.kappa = [];
end
% weight for causal occlusion constraint
if ~isempty(past.layers);
  problem.init_layers = past.t0.layers(:);
end

tic;
layers = pd_wrapper(problem);
layers = reshape(double(layers), imsize);
toc;

%------------------------------------------------------------------------
% TODO: filter the layers
%------------------------------------------------------------------------
layers = postfilter_layers(layers, i1_bflt, 0.5);
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
% objects
%------------------------------------------------------------------------
[object_map_snap, last_obj_id] = layers_to_detachable_objects(layers);
[objects_out, object_map_out] = tao_update_objects( ...
  objects, object_map, object_map_snap, uvb_bflt, uvf_bflt);

objects = objects_out;
object_map = object_map_out;

%------------------------------------------------------------------------
% visualization & output
%------------------------------------------------------------------------
tao_visual(layers, ...
  constraints_causal_b, constraints_causal_f, constraint_weights_old_b, ...
  constraint_weights_old_f, constraints_now_b, constraints_now_f, ...
  constraint_weights_now_nodiv_b, constraint_weights_now_nodiv_f, ...
  constraint_weights_now_b, constraint_weights_now_f, constraints, ...
  constraint_weights, ...
  occb, occf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob, ...
  occb_rev_prob, occf_rev_prob, ...
  unary_constraints_map, max(0.0, 1.0 - weights_now), ...
  max(0.0, 1.0 - weights), imsize, uvb, uvf, ...
  uvb_bflt, uvf_bflt, uvb_cbf_bflt, uvf_cbf_bflt, ...
  uvb_rev, uvf_rev, uvb_cbf, uvf_cbf, ...
  prob_fg, i1, i1_bflt, I1, E, outpath, nameStr, params.versiontype, ...
  seq, k, FINISH, past, params.OCCPROB, problem, boxes, ...
  params.BOX_RAD, params.BOXHELP, params.VIS, params, ...
  objects, object_map, Wx, weights, weights_now);
end
