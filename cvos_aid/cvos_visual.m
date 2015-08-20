%------------------------------------------------------------------------
% cvos_visual
%
% provides all visualizations for our code and stores results
%------------------------------------------------------------------------
function cvos_visual(layers, ...
  constraints_causal_b, constraints_causal_f, constraint_weights_old_b, ...
  constraint_weights_old_f, constraints_now_b, constraints_now_f, ...
  constraint_weights_now_nodiv_b, constraint_weights_now_nodiv_f, ...
  constraint_weights_now_b, constraint_weights_now_f, constraints, ...
  constraint_weights, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob, ...
  weights_now, weights, wx_vis, imsize, uvb_cbf, uvf_cbf, ...
  prob_fg, i1, i1_bflt, I1, outpath, nameStr, versiontype, seq, k, ...
  past, problem, boxes, BOX_RAD, BOXHELP, VIS, params, ...
  object_map, kappa_fg, kappa_box)

rows = imsize(1);
cols = imsize(2);
n = rows * cols;

SAVEFIGBOX = true;

%------------------------------------------------------------------------
% visualization 
%------------------------------------------------------------------------
% draws the layer legend
lay = reshape(layers, imsize);
layermap = lay;
layermap(1,1,:)=0;
layermap(1,2,:)=3;
layermap(1:10,1:10)=0;
layermap(1:10,11:20)=1;
layermap(1:10,21:30)=2;
layermap(1:10,31:40)=3;
layermap = min(3, layermap); % cuts off anything higher

layerimg = 0.4 * sc(i1, 'gray') + 0.6 * sc(layermap, 'jet');
% most frames have 6 or less objects of interest
objectimg = 0.4 * sc(i1, 'gray') + 0.6 * sc(object_map, 'jet', [0,4]);

%------------------------------------------------------------------
% prob_fg
%------------------------------------------------------------------
prob_fg_visual = prob_fg;
prob_fg_visual_max = max(vec(prob_fg_visual));
prev_layer_img = 0.4 * sc(i1, 'gray') + 0.6 * ...
  sc(prob_fg_visual, 'jet', [0.0, 5.0]);

%------------------------------------------------------------------
% visualization and output
%------------------------------------------------------------------
% aggregating constraints
[constraints_causal, constraint_weights_causal] = aggregate_pairs_fast( ...
  double([constraints_causal_b; constraints_causal_f]), ...
  double([constraint_weights_old_b; constraint_weights_old_f]));

constraints_now = aggregate_pairs_fast( ...
  double([constraints_now_b; constraints_now_f]), ...
  double([constraint_weights_now_nodiv_b; ...
  constraint_weights_now_nodiv_f]));

[~, constraint_weights_now] = aggregate_pairs_fast( ...
  double([constraints_now_b; constraints_now_f]), ...
  double([constraint_weights_now_b; constraint_weights_now_f]));

% creating constraints / weights images
past_weights = past.t0.weights;
if isempty(past_weights); past_weights = weights; end;
vis_weights_old = vis_weights_occlusions(max(0.0, 1.0 - past_weights), ...
  constraints_causal, imsize, constraint_weights_causal);
vis_weights_now = vis_weights_occlusions(weights_now, constraints_now, ...
  imsize, constraint_weights_now);
vis_weights_final = vis_weights_occlusions(weights, constraints, ...
  imsize, constraint_weights);

vis_weight_img = reshape(weights(1:n) + weights((n+1):(2*n)), imsize) / 2;
uvb_cbf_img = im2double(flowToColor(uvb_cbf));
uvf_cbf_img = im2double(flowToColor(uvf_cbf));

kappa = problem.kappa;
if isempty(kappa); kappa = zeros(imsize); end;  
kappa_img = sc(reshape(clip(kappa, 0, 1), imsize), 'jet', [0, params.PROB_FG]);
kappa_img = 0.6 * kappa_img + 0.4 * repmat(i1, [1,1,3]);

if (isempty(kappa_box) || (numel(kappa_box) < n)); kappa_box = zeros(imsize); end;
kappa_box_img = sc(reshape(clip(kappa_box, 0, 1), imsize), 'jet', [0, params.PROB_FG]);
kappa_box_img = 0.6 * kappa_box_img + 0.4 * repmat(i1, [1,1,3]);

if isempty(kappa_fg); kappa_fg = zeros(imsize); end;  
kappa_fg_img = sc(reshape(clip(kappa_fg, 0, 1), imsize), 'jet', [0, params.PROB_FG]);
kappa_fg_img = 0.6 * kappa_fg_img + 0.4 * repmat(i1, [1,1,3]);

% drawing the images
if VIS < 500;
  fig(150); clf;
  %------------------------------------------------------------------
  % TOP ROW
  %------------------------------------------------------------------
  vl_tightsubplot(3,4,1);
  imagesc(vis_weights_old); axt; notick;
  th = text(cols/2 - cols/3,rows/10,sprintf('constraints: casual: %0.2f', ...
    max(constraint_weights_causal)));
  set(th, 'Color', 'blue');

  vl_tightsubplot(3,4,2);
  imagesc(vis_weights_now); axt; notick; title(k);
  th = text(cols/2 - cols/5,rows/10,sprintf('current: %0.2f', ...
    max(constraint_weights_now)));
  set(th, 'Color', 'blue');

  vl_tightsubplot(3,4,3);
  imagesc(vis_weights_final); axt; notick;
  th = text(cols/2 - cols/5,rows/10,sprintf('combined: %0.2f', ...
    max(constraint_weights)));
  set(th, 'Color', 'blue');

  vl_tightsubplot(3,4,4);
  imagesc(vis_weight_img, [0, max(vis_weight_img(:))]); axt; notick;
  th = text(cols/2 - cols/5,rows/10,sprintf('weights (jet): %0.2f', ...
    max(problem.Wx(:))));
  set(th, 'Color', 'white');

  %------------------------------------------------------------------
  % MID ROW
  %------------------------------------------------------------------
  vl_tightsubplot(3,4,5);
  imagesc(i1_bflt); axt; notick; title(k);
  
  vl_tightsubplot(3,4,6);
  th = text(cols/2 - cols/5,rows/10,sprintf('weights (jet): %0.2f', ...
    max(vis_weight_img(:))));
  set(th, 'Color', 'white');
  imagesc(uvb_cbf_img); axt; notick;
  vl_tightsubplot(3,4,7);
  th = text(cols/2 - cols/5,rows/10,sprintf('weights (jet): %0.2f', ...
    max(vis_weight_img(:))));
  set(th, 'Color', 'white');
  imagesc(uvf_cbf_img); axt; notick;

  vl_tightsubplot(3,4,8);
  imagesc(prev_layer_img); axt; notick;
  th = text(cols/2 - cols/5,rows/10,sprintf('fg prior: %0.2f', ...
    prob_fg_visual_max));
  set(th, 'Color', 'white');

  %------------------------------------------------------------------
  % BOT ROW
  %------------------------------------------------------------------
  vl_tightsubplot(3,4,9);
  imagesc(layerimg); axt; notick; % output layers
  th = text(cols/2 - cols/5, rows/10, sprintf('layers: %0.2f', ...
    max(layers(:))));
  set(th, 'Color', 'white');
  vl_tightsubplot(3,4,10);
  imagesc(occb_cbf_prob); axt; notick;
  vl_tightsubplot(3,4,11);
  imagesc(occf_cbf_prob); axt; notick;

  %------------------------------------------------------------------
  % Other images
  %------------------------------------------------------------------
  % kappa (with newly obtained boxes)
  fig(150);
  vl_tightsubplot(3,4,12);
  tt = sprintf('kappa: %0.2f', max(vec(kappa)));
  if BOXHELP && ~isempty(boxes);
    box_xy = bbox_to_box_center(boxes);
    box_conf = cat(1, boxes.conf);
    [~, box_img] = draw_boxes(layers, box_xy, BOX_RAD, box_conf);
    b = find(box_img > 0);
    kappa_img(b) = 0.5 * kappa_img(b) + 0.5;
    kappa_img(b + n) = 0.5 * kappa_img(b + n) + 0.5;
    kappa_img(b + 2*n) = 0.5 * kappa_img(b + 2*n) + 0.5;
    tt = strrep(tt, 'kappa', 'kappa (+boxes)');
  end
  imagesc(kappa_img); axt; notick;
  th = text(cols/2 - cols/5, rows/10, tt);
  set(th, 'Color', 'white');

  drawnow;
end

if SAVEFIGBOX;
 %--------------- boxes figure ---------------------
  % boxes + fg prior
  % kappa (with old boxes for this segmentation)
  fgkappa = kappa_fg_img;
  bkappa = kappa_box_img;
  fgbkappa = sc(reshape(clip(kappa, 0, 1), imsize), 'jet', [0, params.PROB_FG]);
  fgbkappa = 0.6 * fgbkappa + 0.4 * repmat(i1, [1,1,3]);
  tt = sprintf('kappa: %0.2f', max(vec(kappa)));
  if BOXHELP && ~isempty(past.boxes);
    pbox_xy = bbox_to_box_center(past.boxes);
    pbox_conf = cat(1, past.boxes.conf);
    [~, pbox_img] = draw_boxes(layers, pbox_xy, BOX_RAD, pbox_conf);
    pb = find(pbox_img > 0);

    fgkappa(pb) = 0.5 * fgkappa(pb) + 0.5;
    fgkappa(pb + n) = 0.5 * fgkappa(pb + n) + 0.5;
    fgkappa(pb + 2*n) = 0.5 * fgkappa(pb + 2*n) + 0.5;
    
    bkappa(pb) = 0.5 * bkappa(pb) + 0.5;
    bkappa(pb + n) = 0.5 * bkappa(pb + n) + 0.5;
    bkappa(pb + 2*n) = 0.5 * bkappa(pb + 2*n) + 0.5;
    
    fgbkappa(pb) = 0.5 * fgbkappa(pb) + 0.5;
    fgbkappa(pb + n) = 0.5 * fgbkappa(pb + n) + 0.5;
    fgbkappa(pb + 2*n) = 0.5 * fgbkappa(pb + 2*n) + 0.5;
  end
    
  if VIS < 400;
    fig(88); clf; imagesc([fgbkappa; kappa_box_img; kappa_fg_img]); 
    axt; notick; title(tt);
  end
  name = sprintf(nameStr, seq, k);
 
  % % uncomment to save a struct with a lot more info than just the result
  % A = v2struct(fgkappa, bkappa, fgbkappa, kappa_img, i1, problem, ...
  %   object_map, layermap, lay, layerimg, kappa_fg_img, kappa_box_img, ...
  %   prob_fg);
  % Amat = fullfile(outpath, [name, '_fgs_', versiontype, '.mat']);
  % save(Amat, '-struct', 'A');
  
  % % uncomment to print fg prior after applying local shape classifiers
  % aafgb = fullfile(outpath, [name, '_fgb_', versiontype, '.png']);
  % imwrite(uint8(255 * kappa_img), aafgb, 'png');

  % % uncomment to print fg prior
  % aafg = fullfile(outpath, [name, '_fg_', versiontype, '.png']);
  % imwrite(uint8(255 * kappa_fg_img), aafg, 'png');
  
  % % uncomment to save weights images
  % wxn = min(clip(weights_now, 0.0, 1.0), [], 3);
  % wxf = min(clip(wx_vis.Wx, 0.0, 1.0), [], 3);
  % aawx = fullfile(outpath, [name, '_wxn_old_', versiontype, '.png']);
  % imwrite(uint8(255 * sc(wxn, 'gray')), aawx, 'png'); % old weights now
  % aawx = fullfile(outpath, [name, '_wxf_old_', versiontype, '.png']);
  % imwrite(uint8(255 * sc(wxf, 'gray')), aawx, 'png'); % old weights
  % wxn = max(wx_vis.weights_now, [], 3);
  % wxf = max(wx_vis.weights, [], 3);
  % wxu = max(wx_vis.unity, [], 3);
  % pwx = min(wx_vis.Wx, [], 3);
  % 
  % % fig grayscale
  % aawx = fullfile(outpath, [name, '_wxn_', versiontype, '.png']);
  % imwrite(uint8(255*sc(clip(1.0 - wxn, 0, 1), 'gray')), aawx, 'png');
  % aawx = fullfile(outpath, [name, '_wxf_', versiontype, '.png']);
  % imwrite(uint8(255*sc(1.0 - wxf, 'gray')), aawx, 'png');
  % aawx = fullfile(outpath, [name, '_wxu_', versiontype, '.png']);
  % imwrite(uint8(255*sc(wxu, 'gray')), aawx, 'png');
  % aawx = fullfile(outpath, [name, '_pwx_', versiontype, '.png']);
  % imwrite(uint8(255*sc(pwx, 'gray')), aawx, 'png');
  %
  % % fig jet
  % aawx = fullfile(outpath, [name, '_wxn_c_', versiontype, '.png']);
  % imwrite(uint8(255*sc(clip(wxn, 0, 1), 'jet')), aawx, 'png');
  % aawx = fullfile(outpath, [name, '_wxf_c_', versiontype, '.png']);
  % imwrite(uint8(255*sc(wxf, 'jet')), aawx, 'png');
  % aawx = fullfile(outpath, [name, '_wxu_c_', versiontype, '.png']);
  % imwrite(uint8(255*sc(wxu, 'jet')), aawx, 'png');
  % aawx = fullfile(outpath, [name, '_pwx_c_', versiontype, '.png']);
  % imwrite(uint8(255*sc(pwx, 'jet')), aawx, 'png');
end

%--------------------------------------------------------------------
% Saving results / printing
%--------------------------------------------------------------------
name = sprintf(nameStr, seq, k);
outResFile = fullfile(outpath, name);

% printing result image
outImgFile = fullfile(outpath, [name, '_img_', versiontype, '.png']);
imwrite(imresize(layerimg, 0.5), outImgFile, 'png');

% object image
outImgFile = fullfile(outpath, [name, '_obj_', versiontype, '.png']);
imwrite(imresize(objectimg, 0.5), outImgFile, 'png');

% printing cues visualization (occluders + occlusions on original image)
vis_cues_img = vis_cues(I1, constraints, min(1.0, constraint_weights));
outCuesFile = fullfile(outpath, [name, '_cue_', versiontype, '.png']);
imwrite(imresize(vis_cues_img, 0.5), outCuesFile, 'png');

%------------------------------------------------------------------
% results saving
%------------------------------------------------------------------
save([outResFile, '_lay_', versiontype, '.mat'], 'lay', 'object_map');
end
