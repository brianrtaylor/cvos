%-----------------------------------------------------------------------------
% cvos(DATA, PKG, ADD)
%
% causal video object segmentation from persistence of occlusions
%
% @param: DATA : name of the sequence to run
% @param: PKG : struct of paths, settings, and parameters
% @param: ADD : struct with data for processing the first and last frames
%-----------------------------------------------------------------------------
function [out, etc] = cvos(DATA, PKG, ADD)
%-----------------------------------------------------------------------------
out = []; etc = [];
if exist('ADD','var');
  %---------------------------------------------------------------------------
  % special case first/last frame setup
  %---------------------------------------------------------------------------

  % obtains params, seq, flow_path, img_path, outpath, files, flow_files,
  % edge_model, past, boxes, dx_inds, dy_inds, Dx, Dy, problem, nameStr, 
  % object_mean_uvf_map, T, FXF
  kk = PKG;
  v2struct(ADD);
  
  fprintf('=========================== time: %d =============\n', kk)
  BEGIN = kk-1;
  LAST = BEGIN > FINISH;
  if LAST; BEGIN = 0; end;
  
else
  %---------------------------------------------------------------------------
  % typical startup
  %---------------------------------------------------------------------------
  BEGIN = 1;
  %---------------------------------------------------------------------------
  % data path
  %---------------------------------------------------------------------------
  [seq, flow_path, img_path, ~, extension, flowtype] = dataPaths(DATA);
  files = dir([img_path '/*.' extension]);
  flow_files = dir([flow_path '/*.mat']);
  
  T = length(files);
  I1 = imread(fullfile(img_path, files(BEGIN).name));
  I2 = imread(fullfile(img_path, files(BEGIN+1).name));
  if (size(I1, 3) == 1), I1 = repmat(I1, [1 1 3]); end;
  if (size(I2, 3) == 1), I2 = repmat(I2, [1 1 3]); end;
  [rows, cols, ~] = size(I1);
  imsize = [rows, cols];
  uvsize = [rows, cols, 2];
  [XX, YY] = meshgrid(1:cols, 1:rows);
  i1 = im2double(rgb2gray(I1));
  i2 = im2double(rgb2gray(I2));
  
  %---------------------------------------------------------------------------
  % algorithm settings
  %---------------------------------------------------------------------------
  params = cvos_params_default();
  edge_model = load('modelFinal.mat');
  
  %---------------------------------------------------------------------------
  % overwrite any variables used in here before running
  %---------------------------------------------------------------------------
  if exist('PKG','var') && isstruct(PKG); 
    params = structmerge(params, PKG); 
  end
  %---------------------------------------------------------------------------
  if ~isfield(params, 'TEST'); params.TEST = false; end;
  if params.DO_FORBACKCAUSAL;
    files = cat(1, files, files((T-1):-1:1));
    T = length(files);
  end
  FXF = strfind(params.model, 'fxf');
  FXF = ~isempty(FXF) && FXF == 1;
  
  %--------------------------------------------------------------------------
  % output files
  %--------------------------------------------------------------------------
  outpath = [params.outpath, '/', seq];
  if ~exist(outpath,'dir'), mkdir(outpath); end;
  out_fname = fullfile(outpath, sprintf('%s_cvosresults.mat', seq));
  out_working_fname = [out_fname, '.running'];
  unix(sprintf('touch %s', out_working_fname));
  
  %--------------------------------------------------------------------------
  % Already done, so don't recompute
  %--------------------------------------------------------------------------
  if ~exist('nameStr', 'var'); nameStr = '%s_%06d'; end;
  out_lay_file = fullfile(outpath, sprintf('%s_lay.mat', seq));
  out_obj_file = fullfile(outpath, sprintf('%s_obj.mat', seq));
  out_fb_lay_file = fullfile(outpath, sprintf('%s_fb_lay.mat', seq));
  out_fb_obj_file = fullfile(outpath, sprintf('%s_fb_obj.mat', seq));
  
  if ((params.DO_FORBACKCAUSAL ...
    && exist(out_fname, 'file') ...
    && exist(out_lay_file, 'file') ...
    && exist(out_obj_file, 'file') ...
    && exist(out_fb_lay_file, 'file') ...
    && exist(out_fb_obj_file, 'file')) ...
    || (~params.DO_FORBACKCAUSAL ...
    && exist(out_fname, 'file') ...
    && exist(out_lay_file, 'file')));
    fprintf('%s: "%s" already finished\n', mfilename, seq);
    fprintf('results files are in: %s\n', outpath);
    return;
  end
  %-----------------------------------------------------------------------
  % first frame init, for storying "history"
  %-----------------------------------------------------------------------
  past = struct;
  past.uvf_cbf = [];
  past.uvf = [];
  past.uvf_rev = [];
  
  past.uvb = [];
  past.uvb_rev = [];
  past.uvf = [];
  past.uvf_rev = [];
  
  past.occf_cbf = [];
  past.occf_rev = [];
  
  past.constraints_b = [];
  past.constraint_weights_b = [];
  past.constraints_causal_b = [];
  past.constraint_weights_causal_b = [];
  past.constraints_f = [];
  past.constraint_weights_f = [];
  past.constraints_causal_f = [];
  past.constraint_weights_causal_f = [];
  
  past.weights = [];
  past.constraints = [];
  past.constraint_weights = [];
  past.constraint_ages = [];
  past.prob_fg = [];
  past.xi_b = [];
  past.xi_f = [];
  past.layers = [];
  past.prob_fg = [];
  past.prob_fg_layer = [];
  
  past.unity = [];
  past.unity_layer = [];
  
  % things warped to current domain in past, are called past.t0.
  past.t0.layers = [];
  past.t0.weights = [];
  
  %-----------------------------------------------------------------------
  % adding to params struct
  %-----------------------------------------------------------------------
  GMMPKG = struct;
  GMMPKG.cons_perturb = params.CONS_PERTURB;
  params.GMMPKG = GMMPKG;
  params.edge_model = edge_model;
  
  %-----------------------------------------------------------------------
  % problem setup
  %-----------------------------------------------------------------------
  problem = struct;
  problem.sigma_c = 0.90/sqrt(8.0);
  problem.sigma_y = 0.90/sqrt(8.0);
  problem.max_iterations = 8000;
  problem.verbosity = 5000;
  problem.fx_tolerance = 0;
  problem.dx_tolerance = 1e-9;
  problem.SOLVE_PIXELWISE = 1;
  problem.layer_upper_bound = 3;
  problem.imsize = imsize;
  problem.nnodes = prod(imsize(1:2));
  problem.nedges = 2*prod(imsize(1:2));
  problem.tau1 = params.TAU1 * ones(problem.nnodes, 1);
  
  %-----------------------------------------------------------------------
  % boxes
  %-----------------------------------------------------------------------
  % useful values
  [Dx, Dy, dx_inds, dy_inds] = make_difference_operator(imsize);
  
  % using boxes for segmentation
  boxes = []; I1lab = []; objects = [];
  object_map = zeros(imsize);
  
  processed_prior_frame = 1;
  checkpoint_timer = params.CHECKPOINT;
  
  %-----------------------------------------------------------------------
  % load from existing work if a causal model is in use
  %-----------------------------------------------------------------------
  if params.CAUSAL && ~FXF;
    checkpointFileName = sprintf('%s_*.mat', out_fname);
    checkpointFiles = dir(checkpointFileName);
    
    if ~isempty(checkpointFiles);
      lastFile = fullfile(outpath, checkpointFiles(end).name);
      
      try
        checkpoint_data = load(lastFile);
        
        % important variables
        past = checkpoint_data.past;
        boxes = checkpoint_data.boxes;
        BEGIN = checkpoint_data.k;
        objects = checkpoint_data.objects;
        object_map = checkpoint_data.object_map;
        
        I1 = checkpoint_data.I1;
        I2 = checkpoint_data.I2;
        i1 = im2double(rgb2gray(I1));
        i2 = im2double(rgb2gray(I2));
        I1_bflt = recursive_bf_mex(I1, 0.01, 0.1, 1, 5);
        I1lab = vl_xyz2lab(vl_rgb2xyz(uint8(I1_bflt)));
      catch e
        fprintf('checkpoint file corrupted, please remove\n%s\n', lastFile);
        display(e);
        return
      end
    end
  end
end

%=============================================================================
% run the sequence
%=============================================================================
BEGIN = BEGIN + 1;
FINISH = T - 1;
if BEGIN > FINISH; FINISH = BEGIN; end;
for k = BEGIN:FINISH;
  fprintf('================== time: %d / %d ===============\n', k, FINISH);
  %------------------------------------------------------------------
  % prevent repeated work for fxf model
  %------------------------------------------------------------------  
  if FXF && k < FINISH;    
    outResFile = fullfile(outpath, sprintf(nameStr, seq, k));
    if exist([outResFile, '_lay_', params.versiontype, '.mat'], 'file');
      processed_prior_frame = 0;
      continue;
    end
    
    % load appropriate stuff
    if processed_prior_frame == 0;
      I1 = imread(fullfile(img_path, files(k).name));
      I2 = imread(fullfile(img_path, files(k+1).name));
      if (size(I1,3)==1), I1 = repmat(I1,[1 1 3]); end;
      if (size(I2,3)==1), I2 = repmat(I2,[1 1 3]); end;
      i1 = im2double(rgb2gray(I1));
      i2 = im2double(rgb2gray(I2));
    end

    % since we have prior frame setup correctly
    processed_prior_frame = 1;
  end

  %------------------------------------------------------------------
  % read 1 new image + pre-processing
  %------------------------------------------------------------------
  a = tic();

  if ~exist('ADD', 'var');
    % transition from t --> t+1
    [I0, I1, i0, i1] = deal(I1, I2, i1, i2);
    I2 = imread(fullfile(img_path, files(k+1).name));
    if (size(I2,3)==1), I2 = repmat(I2,[1 1 3]); end;
    i2 = im2double(rgb2gray(I2));
    I0lab = I1lab;

    % pre-processing on current t
    I1_bflt = recursive_bf_mex(I1, 0.01, 0.1, 1, 5);
    i1_bflt = im2double(I1_bflt);
    I1lab = vl_xyz2lab(vl_rgb2xyz(uint8(I1_bflt)));
  else
    checkpoint_timer = 10000; % big constant so not used
    if k == 1;
      I1 = imread(fullfile(img_path, files(k).name));
      I2 = imread(fullfile(img_path, files(k + 1).name));
      if (size(I1,3)==1), I1 = repmat(I1,[1 1 3]); end;
      if (size(I2,3)==1), I2 = repmat(I2,[1 1 3]); end;

      [rows, cols, ~] = size(I1);
      imsize = [rows,cols]; uvsize = [rows,cols,2];
      i1 = im2double(rgb2gray(I1)); i2 = im2double(rgb2gray(I2));
      I0 = zeros(rows, cols, 3); i0 = zeros(imsize); I0lab = I0;

      % pre-processing on current t
      I1_bflt = recursive_bf_mex(I1, 0.01, 0.1, 1, 5);
      i1_bflt = im2double(I1_bflt);
      I1lab = vl_xyz2lab(vl_rgb2xyz( uint8(i1_bflt*255) ));
    else
      I1 = imread(fullfile(img_path, files(k).name));
      I0 = imread(fullfile(img_path, files(k - 1).name));
      if (size(I1,3)==1), I1 = repmat(I1,[1 1 3]); end;
      if (size(I0,3)==1), I0 = repmat(I0,[1 1 3]); end;

      [rows, cols, ~] = size(I1);
      imsize = [rows,cols]; uvsize = [rows,cols,2];
      i1 = im2double(rgb2gray(I1)); i0 = im2double(rgb2gray(I0));
      I2 = zeros(rows, cols, 3); i2 = zeros(imsize); I2lab = I2; 

      % pre-processing on current t
      I1_bflt = recursive_bf_mex(I1, 0.01, 0.1, 1, 5);
      i1_bflt = im2double(I1_bflt);
      I1lab = vl_xyz2lab(vl_rgb2xyz( uint8(i1_bflt*255) ));
      I0_bflt = recursive_bf_mex(I0, 0.01, 0.1, 1, 5);
      i0_bflt = im2double(I0_bflt);
      I0lab = vl_xyz2lab(vl_rgb2xyz( uint8(i0_bflt*255) ));
    end
  end

  
  E = edgesDetect(I1_bflt, edge_model.model);
  EdgeDollar = 0.75*(1-E);
  
  t_preprocess = toc(a);
  fprintf('C: preprocessing images: %0.3f\n', t_preprocess);

  %------------------------------------------------------------------
  % flow + occlusions + pre-processing
  %------------------------------------------------------------------
  a = tic();
  [uvb, uvb_cbf, uvf, uvf_cbf, uvb_rev, uvf_rev, ...
    ~, occb_cbf, occb_cbf_prob, ~, occf_cbf, occf_cbf_prob, ~, occf_rev] = ...
    cvos_flow_occ(flow_path, flow_files, seq, k, T, ...
    i0, i1, i2, I0, I1, I2, past, params);
  
  t_flowocc_cbf = toc(a);
  fprintf('C: calc flow and occlusions (cbf): %0.3f\n', t_flowocc_cbf);
  a = tic();

  uvf_bflt = recursive_bf_mex(double(uvf), 0.05, 0.004, 0, 10);
  uvb_bflt = recursive_bf_mex(double(uvb), 0.05, 0.004, 0, 10);
  uvf_cbf_bflt = recursive_bf_mex(double(uvf_cbf), 0.05, 0.004, 0, 10);
  uvb_cbf_bflt = recursive_bf_mex(double(uvb_cbf), 0.05, 0.004, 0, 10);
  uvf_rev_bflt = recursive_bf_mex(double(uvf_rev), 0.05, 0.004, 0, 10);
  uvb_rev_bflt = recursive_bf_mex(double(uvb_rev), 0.05, 0.004, 0, 10);
  
  t_flowocc_bf = toc(a);
  fprintf('C: calc flow and occlusions (bf): %0.3f\n', t_flowocc_bf);
  
  %------------------------------------------------------------------
  % layers put into the current frame reference frame
  %------------------------------------------------------------------
  a = tic();
  
  occb_mask = occb_cbf_prob > params.OCCPROBLAYERTHRESH;
  
  if ~isempty(past.layers) && params.CAUSAL;
    fprintf('removing occb mask.\n')
    layers_t0 = round(utils_warp_image(past.layers, uvb_cbf_bflt));
    [past.t0.layers] = remove_occb_mask( layers_t0, occb_mask, 0, 0);
  end
  
  t_layersnow = toc(a);
  fprintf('C: layers warped to the current frame: %0.3f\n', t_layersnow);

  %------------------------------------------------------------------
  % boxes to probability images
  %------------------------------------------------------------------
  a = tic();
  
  if ~isempty(past.layers) && params.CAUSAL;
    fprintf('C: computing local shape classifiers.\n');  
    [boxes, prob_box_fg, ~, ~, ~, weights_box] = cvos_get_prob_from_boxes( ...
      occb_mask, past, I1lab, boxes, params);
  end
  
  t_box_probs = toc(a);
  fprintf('C: compute prob box images: %0.3f\n', t_box_probs);
  
  %------------------------------------------------------------------
  % constraints now
  %------------------------------------------------------------------
  a = tic();
  
  [constraints_now_b, constraint_weights_now_b] ...
    = make_cvos_constraints(occb_cbf_prob, uvb_cbf_bflt, uvb_rev_bflt, ...
    params.OCCPROB, params.MINCONSTRAINTDIST, params.VIS, 16);
  [constraints_now_f, constraint_weights_now_f] ...
    = make_cvos_constraints(occf_cbf_prob, uvf_cbf_bflt, uvf_rev_bflt, ...
    params.OCCPROB, params.MINCONSTRAINTDIST, params.VIS);
  
  t_constraints_now = toc(a);
  fprintf('C: constraints now: %0.3f\n', t_constraints_now);

  %------------------------------------------------------------------
  % weights now
  %------------------------------------------------------------------
  a = tic();
  
  [weights_now, weights_inds] = make_cvos_weights_current_frame( ...
    I1lab, uvb_cbf_bflt, uvf_cbf_bflt, dx_inds, dy_inds, ...
    EdgeDollar, params);
  
  t_weights_now = toc(a);
  fprintf('C: weights now: %0.3f\n', t_weights_now);  
  
  %------------------------------------------------------------------
  % constraints weights now
  %------------------------------------------------------------------
  a = tic();
  if ~isempty(constraints_now_b);
    [constraints_now_b, constraint_weights_now_b, ...
      constraint_weights_now_nodiv_b] = make_cvos_constraint_weights( ...
      constraints_now_b, constraint_weights_now_b, uvf_bflt, uvb_cbf_bflt, ...
      I1lab, weights_now, EdgeDollar, params, 18);
  else
    constraint_weights_now_nodiv_b = [];
  end
  
  if ~isempty(constraints_now_f);
    [constraints_now_f, constraint_weights_now_f, ...
      constraint_weights_now_nodiv_f] = make_cvos_constraint_weights( ...
      constraints_now_f, constraint_weights_now_f, uvb_bflt, uvf_cbf_bflt, ...
      I1lab, weights_now, EdgeDollar, params, 19);
  else
    constraint_weights_now_nodiv_f = [];
  end
  
  t_constraint_weights_now = toc(a);
  fprintf('C: constraint weights now: %0.3f\n', t_constraint_weights_now);
  
  %------------------------------------------------------------------
  % weights from prior frame
  %------------------------------------------------------------------
  a = tic();
  
  if ~isempty(past.weights) && params.CAUSAL && params.WEIGHTHELP;
    [weights, wx_l, wy_l, past] = propagate_cvos_weights( ...
      weights_now, uvb_cbf_bflt, past, occb_mask, ...
      Dx, Dy, k, params);
  else
    weights = weights_now;
  end

  t_weights_propagate = toc(a);
  fprintf('C: weights propagated: %0.3f\n', t_weights_propagate);
  
  %---------------------------------------------------------------------
  % Block that introduces a layer unity prior
  %---------------------------------------------------------------------
  a = tic();
  
  if ~isempty(past.unity) && params.UNITYHELP && params.CAUSAL ...
    && params.PROB_UNITY > 0;
    unity = cvos_unity_prior(past, uvb_bflt, uvb_cbf_bflt, ...
      wx_l, wy_l, occb_mask, k, params);
  else
    unity = zeros(uvsize);
  end
  
  t_unity = toc(a);
  fprintf('C: unity prior: %0.3f\n', t_unity);
  
  %------------------------------------------------------------------
  % constraints from prior frame
  %------------------------------------------------------------------
  a = tic();
  
  if ~exist('prob_box_fg', 'var'); prob_box_fg = []; end;
  
  if ~isempty(past.constraints_b) && params.CAUSAL && params.CONSTRAINTHELP;
    [constraints_b, constraint_weights_b, constraints_causal_b, ...
      constraint_weights_old_b] = propagate_cvos_constraint_weights( ...
      I0lab, I1lab, i1_bflt, uvf_bflt, uvb_cbf_bflt, past, ...
      past.constraints_b, past.constraint_weights_b, past.weights, ...
      weights, past.xi_b, EdgeDollar, constraints_now_b, ...
      constraint_weights_now_b, params, objects, 11);
  else
    constraints_b = constraints_now_b;
    constraint_weights_b = constraint_weights_now_b;
    constraints_causal_b = [];
    constraint_weights_old_b = [];
  end

  if ~isempty(past.constraints_f) && params.CAUSAL && params.CONSTRAINTHELP;
    [constraints_f, constraint_weights_f, constraints_causal_f, ...
      constraint_weights_old_f] = propagate_cvos_constraint_weights( ...
      I0lab, I1lab, i1_bflt, uvb_bflt, uvf_cbf_bflt, past, ...
      past.constraints_f, past.constraint_weights_f, past.weights, ...
      weights, past.xi_f, EdgeDollar, constraints_now_f, ...
      constraint_weights_now_f, params, objects, 12);
  else
    constraints_f = constraints_now_f;
    constraint_weights_f = constraint_weights_now_f;
    constraints_causal_f = [];
    constraint_weights_old_f = [];
  end

  % more pruning based on duplicates
  [constraints, constraint_weights, constraint_inds] = ...
    aggregate_pairs_fast(double([constraints_b; constraints_f]), ...
    double([constraint_weights_b; constraint_weights_f]));

  nconstraints_b = size(constraints_b, 1);
  nconstraints_f = size(constraints_f, 1);
  constraint_inds_b = constraint_inds(1:nconstraints_b);
  constraint_inds_f = constraint_inds( ...
    (nconstraints_b+1):(nconstraints_b+nconstraints_f));
  
  t_constraints_propagate = toc(a);
  fprintf('C: constraints propagate: %0.3f\n', t_constraints_propagate);

  %---------------------------------------------------------------------
  % Block that introduces a ``layer-driven foreground prior''
  %---------------------------------------------------------------------    
  a = tic();
  
  if ~isempty(past.prob_fg) && params.CAUSAL && params.PROB_FG > 0;
    prob_fg = cvos_fg_prior(past, uvb_cbf_bflt, occb_mask, params); 
  else
    prob_fg = zeros(imsize);
  end
     
  t_fgprior = toc(a);
  fprintf('C: foreground prior: %0.3f\n', t_fgprior);

  %------------------------------------------------------------------
  % solve problem
  %------------------------------------------------------------------
  a = tic();

  % creates fake constraints if none exist to keep running
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
      ima = [-weights_box(:, :, 1), -weights_box(:, :, 2), ...
        max(-weights_box, [], 3)];
      imb = [weights_b4(:, :, 1), weights_b4(:, :, 2), ...
        max(weights_b4, [], 3)];
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
  
  % visualization
  if params.VIS < 160;
    probwx_map = min(Wx, [], 3);
    probwx_img = sc(clip(probwx_map, 0, 5), 'jet');
    probwx_img = 0.6 * probwx_img + 0.4 * repmat(i1, [1,1,3]);
    fig(6); clf; imagesc(probwx_img);
    title(sprintf('problem.Wx (t = %d): %0.6f', k, max(Wx(:))));
    fig(7); clf; imagesc(probwx_map);
    title(sprintf('problem.Wx (t = %d): %0.6f', k, max(Wx(:))));
  end
  
  if params.VIS < 130; % visualization of weights combination here    
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
  fg_kappa = zeros(imsize);
  box_kappa = zeros(imsize);
  if params.CAUSAL && ~isempty(past.layers) ...
      && ((params.BOXHELP > 0) || (params.PROB_FG > 0));
    problem.USE_TEMPORAL_PENALTY = 1;
    
    if params.CAUSAL && params.BOXHELP && ~isempty(boxes);
      box_kappa = params.PROB_BOX_FG * prob_box_fg;
    else
      box_kappa = 0;
    end
    
    fg_kappa = params.PROB_FG * prob_fg(:);
    problem.kappa = double(max(0.0, fg_kappa(:) + box_kappa(:)));
    
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
  if ~isempty(past.layers) && params.CAUSAL;
    problem.init_layers = past.t0.layers(:);
  end
  
  t_problem_setup = toc(a);
  fprintf('C: setting up the problem: %0.3f\n', t_problem_setup);
  
  a = tic();
  layers = pd_wrapper(problem);
  layers = reshape(double(layers), imsize);
  t_problem_solve = toc(a);
  fprintf('C: solving the problem: %0.3f\n', t_problem_solve);
  
  %------------------------------------------------------------------------
  % filter the layers
  %------------------------------------------------------------------------
  a = tic();  
  layers = postfilter_layers(layers, double(i1_bflt), 0.5);
  t_postfilter = toc(a);
  fprintf('C: postfiltering: %0.3f\n', t_postfilter);
  
  %------------------------------------------------------------------------

  %------------------------------------------------------------------------
  % data to be used in the next iteration
  %------------------------------------------------------------------------
  a = tic();
  
  past.layers               = round(layers);
  past.weights              = weights;
  past.unity                = unity;
  
  past.constraints          = constraints;
  past.constraint_weights   = constraint_weights;
  past.constraints_b        = constraints_b;
  past.constraint_weights_b = constraint_weights_b;
  past.constraints_f        = constraints_f;
  past.constraint_weights_f = constraint_weights_f;

  % flows and such
  past.uvf            = uvf;
  past.uvf_rev        = uvf_rev;
  past.occf_cbf       = occf_cbf;
  past.occf_rev       = occf_rev;
  
  past.uvf_cbf_bflt   = uvf_cbf_bflt;

  mag_uvf = sqrt(sum(uvf_cbf_bflt .* uvf_cbf_bflt, 3));
  mag_uvf_non_trivial = mag_uvf(mag_uvf > params.TINY_UV);
  past.mag_uvf_avg = mean(mag_uvf_non_trivial(:));
  past.w_warp_uvf_denom = max(1.0, past.mag_uvf_avg / 2.0);
  past.w_warp_uvf   = exp( -sqrt(sum(uvf_cbf_bflt .^ 2, 3)) / max( ...
    sqrt(params.WARP_SAFESPEEDSQUARED), past.w_warp_uvf_denom));
  
  %------------------------------------------------------------------------
  % foreground map to be used in the next iteration
  %------------------------------------------------------------------------
  past.prob_fg_layer = utils_get_foreground_map(layers, 0.8);

  if params.RUNNING_FG_PRIOR && params.CAUSAL;
    past.prob_fg = prob_fg;
  elseif params.CAUSAL;
    past.prob_fg = past.prob_fg_layer;
  end

  %------------------------------------------------------------------------
  % extract \xi
  %------------------------------------------------------------------------
  nconstraints = size(constraints, 1);
  if nconstraints > 0;
    Docc = sparse([(1:nconstraints)'; (1:nconstraints)'],...
      [constraints(:,1); constraints(:,2)],...
      [+ones(nconstraints,1); -ones(nconstraints,1)],...
      nconstraints, problem.nnodes); %size of matrix
    past.xi = round(clip(1 - Docc * layers(:), 0, 2));   
    past.xi_b = past.xi(constraint_inds_b);
    past.xi_f = past.xi(constraint_inds_f);
  else
    past.xi = [];
    past.xi_b = [];
    past.xi_f = [];
  end
  
  t_past = toc(a);
  fprintf('C: setting up past for next frame: %0.3f\n', t_past);
  
  
  %------------------------------------------------------------------------
  % learn models for boxes from segmentation (to help on next frame)
  %------------------------------------------------------------------------  
  a = tic();
  
  if params.CAUSAL && params.BOXHELP;  
    past.boxes = boxes;
    boxes = cvos_update_box_models( ...
      boxes, layers, I1lab, Dx, Dy, weights, params);
  end
  
  t_box_models = toc(a);
  fprintf('C: updating box models: %0.3f\n', t_box_models);

  %------------------------------------------------------------------------
  % objects 
  %------------------------------------------------------------------------
  a = tic();
  
  [object_map_snap] = layers_to_detachable_objects(layers);
  [objects, object_map] = cvos_update_objects( ...
    objects, object_map, object_map_snap, uvb_bflt, uvf_bflt);
  
  t_objects_update = toc(a);
  fprintf('C: updating objects: %0.3f\n', t_objects_update);
  
  %------------------------------------------------------------------------
  % objects top pre-processing for mean flow images
  %------------------------------------------------------------------------
  a = tic();
  
  if ~isempty(past.layers) && params.CAUSAL; 
    object_mean_uvf_map = cvos_get_object_flow_maps(objects, object_map);
    
    uvdiff = uvf_cbf_bflt - object_mean_uvf_map;
    fg_uvf_diff = sqrt(sum(uvdiff .^ 2, 3));
    
    past.w_fg_uvf = exp( -fg_uvf_diff / past.w_warp_uvf_denom);
    past.w_fg_uvf(isnan(past.w_fg_uvf)) = 1.0;   
  else
    object_mean_uvf_map = zeros(uvsize);
    past.fg_uvf_diff = zeros(imsize);
    past.w_fg_uvf = ones(imsize);
  end
  
  t_object_uvf = toc(a);
  fprintf('C: updating objects flow maps: %0.3f\n', t_object_uvf);

  %------------------------------------------------------------------------
  % visualization 
  %------------------------------------------------------------------------
  a = tic();
  wx_vis = v2struct(weights_now, weights_cut, ...
    weights_unity, unity, weights, Wx);
  
  cvos_visual(layers, ...
    constraints_causal_b, constraints_causal_f, constraint_weights_old_b, ...
    constraint_weights_old_f, constraints_now_b, constraints_now_f, ...
    constraint_weights_now_nodiv_b, constraint_weights_now_nodiv_f, ...
    constraint_weights_now_b, constraint_weights_now_f, constraints, ...
    constraint_weights, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob, ...
    max(0.0, 1.0 - weights_now), ...
    max(0.0, 1.0 - weights), wx_vis, imsize, uvb_cbf, uvf_cbf, ...
    prob_fg, i1, i1_bflt, I1, outpath, nameStr, params.versiontype, ...
    seq, k, past, problem, boxes, ...
    params.BOX_RAD, params.BOXHELP, params.VIS, params, ...
    object_map, fg_kappa, box_kappa);
  
  t_vis = toc(a);
  fprintf('C: visualization: %0.3f\n', t_vis);
  
  %------------------------------------------------------------------------
  % prevents repeated computation
  %------------------------------------------------------------------------
  if params.CAUSAL && ~FXF && checkpoint_timer <= 0;
    checkpointFileName = sprintf('%s_%06d.mat', out_fname, k);
    save(checkpointFileName, 'k', 'past', 'boxes', ...
      'object_map', 'objects', 'object_mean_uvf_map', 'I1', 'I2');
    checkpoint_timer = PKG.CHECKPOINT;
  end
  checkpoint_timer = checkpoint_timer - 1;
  
  %------------------------------------------------------------------------
  % show timing
  %------------------------------------------------------------------------
  if params.VIS < 1000;
    rest_time = t_preprocess + t_layersnow + t_box_probs ...
      + t_constraints_now + t_weights_now + t_constraint_weights_now ...
      + t_weights_propagate + t_unity + t_constraints_propagate ...
      + t_fgprior + t_problem_setup ...
      + t_postfilter + t_past + t_box_models + t_objects_update ...
      + t_object_uvf + t_vis + t_flowocc_cbf + t_flowocc_bf;
    total_time = rest_time + t_problem_solve;
    fprintf('--------- timing info ---------\n');
    fprintf('total time:              %02.3f\n', total_time);
    fprintf('solver:                  %02.3f\n', t_problem_solve); 
    fprintf('rest:                    %02.3f\n', rest_time);
    fprintf('--------- ----------- ---------\n');
    fprintf('preprocessing images:    %02.3f\n', t_preprocess);
    fprintf('flow & occlusions cbf:   %02.3f\n', t_flowocc_cbf);
    fprintf('flow & occlusions bf:    %02.3f\n', t_flowocc_bf);
    fprintf('layers warped to now:    %02.3f\n', t_layersnow);
    fprintf('calc prob box images:    %02.3f\n', t_box_probs);
    fprintf('constraints now:         %02.3f\n', t_constraints_now);
    fprintf('weights now:             %02.3f\n', t_weights_now);
    fprintf('constraint weight now:   %02.3f\n', t_constraint_weights_now);
    fprintf('weights propagation:     %02.3f\n', t_weights_propagate);
    fprintf('unity prior:             %02.3f\n', t_unity);
    fprintf('constraints weight prop: %02.3f\n', t_constraints_propagate);
    fprintf('foreground prior:        %02.3f\n', t_fgprior);
    fprintf('problem setup:           %02.3f\n', t_problem_setup);
    fprintf('problem solving:         %02.3f\n', t_problem_solve);
    fprintf('postfiltering:           %02.3f\n', t_postfilter);
    fprintf('update past:             %02.3f\n', t_past);
    fprintf('update box models:       %02.3f\n', t_box_models);
    fprintf('update objects:          %02.3f\n', t_objects_update);
    fprintf('update object uvf maps:  %02.3f\n', t_object_uvf);
    fprintf('visualization + saving:  %02.3f\n', t_vis);
    fprintf('----------- ----------- -------\n');    
  end
  
  %--------------------------------------------------------------------------
  % if we should run the first or the last frame
  %--------------------------------------------------------------------------
  if ~params.TEST && ((k == 2) || (k == FINISH)) && ~exist('ADD','var');
    ADD = v2struct(params, ...
      seq, flow_path, img_path, flowtype, outpath, files, flow_files, ...
      edge_model, past, boxes, dx_inds, dy_inds, ...  
      Dx, Dy, problem, nameStr, object_mean_uvf_map, FINISH, T, FXF);

    if (k == 2); % first frame
      cvos(params, k-1, ADD);
    elseif ~params.DO_FORBACKCAUSAL; % end frame 
      cvos(params, k+1, ADD);
    elseif params.DO_FORBACKCAUSAL; % end frame
      cvos(params, k+1, ADD);
    end
    clear ADD; % prevents code from halting on frame 3
  elseif exist('ADD','var');
    fprintf('done\n');
    break;
  end
end

%-----------------------------------------------------------------------------
% save results
%-----------------------------------------------------------------------------
if ~exist('ADD', 'var');
  save(out_fname, 'params');
  unix(sprintf('rm %s', out_working_fname));
  
  try
    if params.DO_FORBACKCAUSAL;
      utils_compress_lay_files(outpath, seq, params.model, 1);
    end
    utils_compress_lay_files(outpath, seq, params.model, 0);
  catch e
    fprintf('%s: error compressing lay, obj file\n', mfilename);
  end
end
end
