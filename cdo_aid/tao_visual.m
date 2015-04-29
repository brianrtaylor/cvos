%------------------------------------------------------------------------
% tao_visual
%
% provides all visualizations for our code and stores results
%------------------------------------------------------------------------
function tao_visual(layers, ...
  constraints_causal_b, constraints_causal_f, constraint_weights_old_b, ...
  constraint_weights_old_f, constraints_now_b, constraints_now_f, ...
  constraint_weights_now_nodiv_b, constraint_weights_now_nodiv_f, ...
  constraint_weights_now_b, constraint_weights_now_f, constraints, ...
  constraint_weights, ...
  occb, occf, occb_cbf, occf_cbf, occb_cbf_prob, occf_cbf_prob, ...
  occb_rev_prob, occf_rev_prob, ...
  unary_constraints_map, weights_now, weights, imsize, ...
  uvb, uvf, uvb_bflt, uvf_bflt, uvb_cbf_bflt, uvf_cbf_bflt, ...
  uvb_rev, uvf_rev, uvb_cbf, uvf_cbf, ...
  prob_fg, i1, i1_bflt, I1, E, outpath, nameStr, versiontype, seq, k, N, ...
  past, OCCPROB, problem, boxes, BOX_RAD, BOXHELP, VIS, params, ...
  objects, object_map, Wx, weights_raw, weights_now_raw);

rows = imsize(1);
cols = imsize(2);
n = rows * cols;

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
% I don't think any image has more than 6 objects of interest
objectimg = 0.4 * sc(i1, 'gray') + 0.6 * sc(object_map, 'jet', [0,4]);

%------------------------------------------------------------------
% prob_fg
%------------------------------------------------------------------
% prob_fg_visual = log(1 + prob_fg);
prob_fg_visual = prob_fg;
prob_fg_visual_max = max(vec(prob_fg_visual));
% prev_layer_img = 0.4 * sc(i1, 'gray') + 0.6 * ...
%   sc(prob_fg_visual, 'jet', [0, max(prob_fg_visual_max, 1.0)]);
prev_layer_img = 0.4 * sc(i1, 'gray') + 0.6 * ...
  sc(prob_fg_visual, 'jet', [0.0, 5.0]);

%------------------------------------------------------------------
% visualization and output
%------------------------------------------------------------------
% aggregating constraints
[constraints_causal, constraint_weights_causal] = aggregate_pairs_fast( ...
  double([constraints_causal_b; constraints_causal_f]), ...
  double([constraint_weights_old_b; constraint_weights_old_f]));

[constraints_now, constraint_weights_now_nodiv] = aggregate_pairs_fast( ...
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
% vis_weights_now_nodiv = vis_weights_occlusions(weights_now, ...
%   constraints_now, imsize, constraint_weights_now_nodiv);
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

kappa_box = problem.kappa_box; 
if (isempty(kappa_box) || (numel(kappa_box) < n)); kappa_box = zeros(imsize); end;
kappa_box_img = sc(reshape(clip(kappa_box, 0, 1), imsize), 'jet', [0, params.PROB_FG]);
kappa_box_img = 0.6 * kappa_box_img + 0.4 * repmat(i1, [1,1,3]);

kappa_fg = problem.kappa_fg; 
if isempty(kappa_fg); kappa_fg = zeros(imsize); end;  
kappa_fg_img = sc(reshape(clip(kappa_fg, 0, 1), imsize), 'jet', [0, params.PROB_FG]);
kappa_fg_img = 0.6 * kappa_fg_img + 0.4 * repmat(i1, [1,1,3]);

% drawing the images
if VIS < 500;
  
% % %   %------------------------------------------------------------------
% % %   % 20141010 fig-constraint-weights: for marathon meeting
% % %   vcn1 = vis_cues(i1, constraints_now, 0.8*ones(size(constraints_now, 1), 1));
% % %   vcn2 = vis_cues(i1, constraints_now, min(1.0, constraint_weights_now_nodiv));
% % %   vcn3 = vis_cues(i1, constraints_now, min(1.0, constraint_weights_now));
% % %   fig(11); clf; imagesc([vcn1; vcn2; vcn3]);
% % %   axt; notick;
% % %   
% % %   CT = 0.25;
% % %   vcn11 = vis_cues(i1, constraints_now, 0.8*ones(size(constraints_now, 1), 1), CT);
% % %   vcn21 = vis_cues(i1, constraints_now, min(1.0, constraint_weights_now_nodiv), CT);
% % %   vcn31 = vis_cues(i1, constraints_now, min(1.0, constraint_weights_now), CT);
% % %   fig(12); clf; imagesc([vcn11; vcn21; vcn31]);
% % %   axt; notick;
% % %   
% % %   outCuesFile3 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesnow_', versiontype, '.png']);
% % %   imwrite(vcn3, outCuesFile3, 'png');
% % %   outCuesFile2 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesnownodiv_', versiontype, '.png']);
% % %   imwrite(vcn2, outCuesFile2, 'png'); 
% % %   outCuesFile1 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesnownoweight_', versiontype, '.png']);
% % %   imwrite(vcn1, outCuesFile1, 'png');
% % %   
% % %   outCuesFile31 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesnowT_', versiontype, '.png']);
% % %   imwrite(vcn31, outCuesFile31, 'png');
% % %   outCuesFile21 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesnownodivT_', versiontype, '.png']);
% % %   imwrite(vcn21, outCuesFile21, 'png'); 
% % %   outCuesFile11 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesnownoweightT_', versiontype, '.png']);
% % %   imwrite(vcn11, outCuesFile11, 'png');
% % %   %------------------------------------------------------------------
  
% % %   %------------------------------------------------------------------
% % %   % 20141010 fig-cbf: for marathon meeting
% % %   
% % %   fig(13); clf;
% % %   flowshow([[uvb; uvb_cbf],[uvf; uvf_cbf]]); axt; notick;  
% % %   
% % %   uvbFile = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %     '_uvb_', versiontype, '.png']);
% % %   imwrite(flowToColor(uvb), uvbFile, 'png');
% % %   uvfFile = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %     '_uvf_', versiontype, '.png']);
% % %   imwrite(flowToColor(uvf), uvfFile, 'png');
% % %   
% % %   uvbCbfFile = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %     '_uvbcbf_', versiontype, '.png']);
% % %   imwrite(flowToColor(uvb_cbf), uvbCbfFile, 'png');
% % %   uvfCbfFile = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %     '_uvfcbf_', versiontype, '.png']);
% % %   imwrite(flowToColor(uvf_cbf), uvfCbfFile, 'png');
% % %   %------------------------------------------------------------------
  
% % %   %------------------------------------------------------------------
% % %   % 20141010 fig-constraints-prob: for marathon meeting
% % %   vcn2 = vis_cues(i1, constraints_now, min(1.0, constraint_weights_now));
% % %   vcn1 = vis_cues(i1, constraints, min(1.0, constraint_weights));
% % %   fig(14); clf; imagesc([vcn2; vcn1]);
% % %   axt; notick;
% % %   
% % %   CT = 0.25;
% % %   vcn21 = vis_cues(i1, constraints_now, min(1.0, constraint_weights_now), CT);
% % %   vcn11 = vis_cues(i1, constraints, min(1.0, constraint_weights), CT);
% % %   fig(15); clf; imagesc([vcn21; vcn11]);
% % %   axt; notick;
% % %   
% % %   outCuesFile2 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesfxf_', versiontype, '.png']);
% % %   imwrite(vcn2, outCuesFile2, 'png'); 
% % %   outCuesFile1 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesintegrated_', versiontype, '.png']);
% % %   imwrite(vcn1, outCuesFile1, 'png');
% % %   
% % %   outCuesFile21 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesfxfT_', versiontype, '.png']);
% % %   imwrite(vcn21, outCuesFile21, 'png'); 
% % %   outCuesFile11 = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %   '_cuesintegratedT_', versiontype, '.png']);
% % %   imwrite(vcn11, outCuesFile11, 'png');
% % %   %------------------------------------------------------------------
% % %   
% % %   %------------------------------------------------------------------
% % %   % 20141010 fig-fg-prior / fig-unity-prior: for marathon meeting
% % %   fg_prior_img = 0.4 * sc(i1, 'gray') + 0.6 * ...
% % %     sc(prob_fg_visual, 'jet', [0.0, 5.0]);
% % %   fig(16); clf; imagesc(fg_prior_img); axt; notick;
% % %   outFgFile = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %     '_fgprior_', versiontype, '.png']);
% % %   imwrite(fg_prior_img, outFgFile, 'png');
% % %   
% % %   Wx = max(0.0, problem.Wx);
% % %   unity_prior = reshape(min(Wx(1:n), Wx((n+1):(2*n))), imsize) / 2;
% % %   unity_prior_raw = min(unity_prior, 5.0);
% % %   unity_img = 0.4 * sc(i1, 'gray') + 0.6 * sc(unity_prior, 'jet', [0.0, 5.0]);
% % %   fig(17); clf; imagesc(unity_img); axt; notick;
% % %   outUnityRawFile = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %     '_unitypriorraw_', versiontype, '.png']);
% % %   outUnityFile = fullfile(outpath, [sprintf(nameStr, seq, k), ...
% % %     '_unitypriorimg_', versiontype, '.png']);
% % %   imwrite(unity_prior_raw, outUnityRawFile, 'png');
% % %   imwrite(unity_img, outUnityFile, 'png');
% % %   %------------------------------------------------------------------

  
  
  fig(150); clf;

  %-----------------------------------
  % TOP ROW
  %-----------------------------------
  vl_tightsubplot(3,4,1);
  imagesc(vis_weights_old); axt; notick;
  th = text(cols/2 - cols/3,rows/10,sprintf('constraints: casual: %0.2f', ...
    max(constraint_weights_causal)));
  set(th, 'Color', 'blue');

  % vl_tightsubplot(3,4,2);
  % imagesc(vis_weights_now_nodiv); axt; notick;
  % th = text(cols/2 - cols/5,rows/10,sprintf('current (nodiv): %0.2f', ...
  %   max(constraint_weights_now_nodiv)));
  % set(th, 'Color', 'blue');

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
%   vis_weight_img = reshape(min(problem.Wx(1:n), problem.Wx((n+1):(2*n))), ...
%     imsize) / 2;
  % vis_weight_img = log(1 + vis_weight_img);
  imagesc(vis_weight_img, [0, max(vis_weight_img(:))]); axt; notick;
  th = text(cols/2 - cols/5,rows/10,sprintf('weights (jet): %0.2f', ...
    max(problem.Wx(:))));
  set(th, 'Color', 'white');

  %-----------------------------------
  % MID ROW
  %-----------------------------------
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

  %-----------------------------------
  % BOT ROW
  %-----------------------------------
  vl_tightsubplot(3,4,9);
  imagesc(layerimg); axt; notick; % output layers
  th = text(cols/2 - cols/5, rows/10, sprintf('layers: %0.2f', ...
    max(layers(:))));
  set(th, 'Color', 'white');
  vl_tightsubplot(3,4,10);
  imagesc(occb_cbf_prob); axt; notick;
  vl_tightsubplot(3,4,11);
  imagesc(occf_cbf_prob); axt; notick;

  % % combined occlusions
  % vl_tightsubplot(3,4,12);
  % imagesc(max(occb_cbf_prob, occf_cbf_prob) > OCCPROB); axt; notick;

  % % unary from occlusion cues (never seemed to help :/, but a good
  % % idea. perhaps we can try again in the future to make it work.
  % % It does detachable objects by unary + smoothness, not pairwise, so
  % % can use the crfs and things)
  % vl_tightsubplot(3,4,12);
  % imagesc(unary_constraints_map, [0, 5]); axt; notick;
  % th = text(cols/2 - cols/5,rows/10,sprintf('range: 0.0 - %0.2f', ...
  %   max(unary_constraints_map(:))));
  % set(th, 'Color', 'white');

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

  drawnow; pause(0.1);
end

SAVEFIGBOX = true;
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
    fig(88); clf; imagesc([fgbkappa; kappa_box_img; kappa_fg_img]); axt; notick; title(tt);
  end
  % th = text(cols/2 - cols/5, rows/10, tt);
  % set(th, 'Color', 'white');
  
  A = v2struct(fgkappa, bkappa, fgbkappa, kappa_img, i1, problem, ...
    object_map, layermap, lay, layerimg, kappa_fg_img, kappa_box_img, ...
    prob_fg);
  
  name = sprintf(nameStr, seq, k);
  Amat = fullfile(outpath, [name, '_fgs_', versiontype, '.mat']);
  save(Amat, '-struct', 'A');
  
  aafgb = fullfile(outpath, [name, '_fgb_', versiontype, '.png']);
  imwrite(uint8(255 * kappa_img), aafgb, 'png');
  aafg = fullfile(outpath, [name, '_fg_', versiontype, '.png']);
  imwrite(uint8(255 * kappa_fg_img), aafg, 'png');
  
  % % % %     aapfgb = fullfile(outpath, [name, '_img_', versiontype, '.pdf']);
  % % % %     imwritepdf(uint8(255 * kappa_img), aapfgb);
  % % % %     aapfg = fullfile(outpath, [name, '_img_', versiontype, '.pdf']);
  % % % %     imwritepdf(uint8(255 * kappa_fg_img), aapfg);
end

%--------------------------------------------------------------------
% Saving results / printing
%--------------------------------------------------------------------
name = sprintf(nameStr, seq, k);
outResFile = fullfile(outpath, name);

% printing display (fig 150)
% gets numbers, but slow / bad for shell
% print(150, '-dpng', [outResFile, '_res_', versiontype, '.png']);
% doesn't get numbers, but faster, and need to change if we change display
% row1 = [vis_weights_old, vis_weights_now, vis_weights_final, ...
%   sc(vis_weight_img, 'jet')];
% row2 = [im2double(I1), uvb_cbf_img, uvf_cbf_img, prev_layer_img];
% row3 = [layerimg, sc(occb_cbf_prob, 'jet'), sc(occf_cbf_prob, 'jet'), kappa_img];
% res_img = [row1; row2; row3];
% imwrite(res_img, [outResFile, '_res_', versiontype, '.png'], 'png');

% printing result image
outImgFile = fullfile(outpath, [name, '_img_', versiontype, '.png']);
imwrite(imresize(layerimg, 0.5), outImgFile, 'png');

% object image
outImgFile = fullfile(outpath, [name, '_obj_', versiontype, '.png']);
imwrite(imresize(objectimg, 0.5), outImgFile, 'png');

% printing cues visualization (occluders + occlusions on original image)
vis_cues_img = vis_cues(I1, constraints, min(1.0, constraint_weights));
% if VIS < 500;
%   fig(151); clf;
%   sc(vis_cues_img); notick;
%   th = text(cols/2 - cols/5,rows/10,sprintf('range: 0.0 - %0.2f', ...
%     max(constraint_weights)));
%   set(th, 'Color', 'blue');
% end
outCuesFile = fullfile(outpath, [name, '_cue_', versiontype, '.png']);
imwrite(imresize(vis_cues_img, 0.5), outCuesFile, 'png');

%------------------------------------------------------------------
% results saving
%------------------------------------------------------------------
%%%% %------------------------------------------------------------------
%%%% % for supplementary material
%%%% %------------------------------------------------------------------
%%%% % layermap      = single(layermap);
%%%% % overlay       = uint8(255*layerimg);
%%%% % vis_weights   = uint8(255*vis_weights_final);
%%%% vis_prob_fg   = single(prob_fg);
%%%% % eb            = single(occb_cbf);
%%%% % ef            = single(occf_cbf);
%%%% % vis_cues_btay = uint8(255*vis_cues_img);
%%%% lay           = uint8(lay);
%%%% object_map    = uint32(object_map);
%%%% constraints   = uint32(constraints);
%%%% constraints_weights  = single(constraint_weights);
%%%% 
%%%% flo.uvf           = single(uvf);
%%%% flo.uvf_cbf       = single(uvf_cbf);
%%%% flo.uvf_bflt      = single(uvf_bflt);
%%%% flo.uvf_cbf_bflt  = single(uvf_cbf_bflt);
%%%% 
%%%% flo.uvb           = single(uvb);
%%%% flo.uvb_cbf       = single(uvb_cbf);
%%%% flo.uvb_bflt      = single(uvb_bflt);
%%%% flo.uvb_cbf_bflt  = single(uvb_cbf_bflt);
%%%% 
%%%% flo.occf          = single(occf);
%%%% flo.occf_cbf      = single(occf_cbf);
%%%% flo.occf_cbf_prob = single(occf_cbf_prob);
%%%% flo.occb          = single(occb);
%%%% flo.occb_cbf      = single(occb_cbf);
%%%% flo.occb_cbf_prob = single(occb_cbf_prob);
%%%%
%%%% save([outResFile, '_lay_', versiontype, '.mat'], ...
%%%%   'lay', 'object_map', 'constraints', 'constraint_weights', ...
%%%%   'vis_prob_fg', 'flo', 'Wx', 'weights_raw', 'weights_now_raw');

%------------------------------------------------------------------
% for original experiments
%------------------------------------------------------------------
layermap      = single(layermap);
overlay       = uint8(255*layerimg);
vis_weights   = uint8(255*vis_weights_final);
vis_prob_fg   = single(prob_fg);
eb            = single(occb_cbf);
ef            = single(occf_cbf);
vis_cues_btay = uint8(255*vis_cues_img);
lay           = uint8(lay);
object_map    = uint32(object_map);
constraints   = uint32(constraints);
constraints_weights  = single(constraint_weights);

save([outResFile, '_lay_', versiontype, '.mat'], 'lay', 'object_map', ...
  'overlay', 'vis_weights', 'vis_prob_fg', 'eb', 'ef', ...
  'vis_cues_btay', 'layermap');
end
