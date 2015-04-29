function fig_cons(seq, id, model, ustr, outpath)
addpath('/home/btay/projects/tao/'); setup;
if ~exist('ustr', 'var'); ustr = 'check'; end;
if ~exist('model', 'var'); model = ''; end;
% if ~exist('outpath', 'var'); outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.6-fullbg-moseg/'; end;
% if ~exist('outpath', 'var'); 
%   outpath = '/plot/btay/projects/detachable/cdov/week18-cdov-r0.5-full/'; 
% end
if ~exist('outpath', 'var'); 
  outpath = '/pad_local/btay/projects/tao/week17/r0p5/'; 
end

nameStr = '%s_%06d';
fname = sprintf(nameStr, seq, id);

FXF = strfind(ustr, 'fxf');

%-------------------------------------------------------------------------
% load data
%-------------------------------------------------------------------------
img_path = '/plot/vasiliy/CVPR15/data/moseg';
f = dir(fullfile(img_path, seq, sprintf('%s_*.png', seq)));
img_file = fullfile(img_path, seq, f(id).name);
f = ls(img_file);

I1 = imread(img_file);
i1 = im2single(I1);
i1g = sc(i1, 'gray');

[rows, cols, ~] = size(I1);
imsize = [rows, cols];

f = dir(fullfile(outpath, seq, sprintf('%s_*%s*_lay_*%s*.mat', ...
  seq, sprintf('%06d', id), model)));
mat_file = fullfile(outpath, seq, f.name);
% mat_file = fullfile(outpath, seq, f(id).name);
A = load(mat_file);

%-------------------------------------------------------------------------
% make images 
%-------------------------------------------------------------------------
% layer image
layermap = A.lay;
% layermap(1,1,:)=0;
% layermap(1,2,:)=3;
% layermap(1:10,1:10)=0;
% layermap(1:10,11:20)=1;
% layermap(1:10,21:30)=2;
% layermap(1:10,31:40)=3;
layermap = min(3, layermap); % cuts off anything higher
lmask = (layermap > 0);
l3mask = repmat(lmask, [1, 1, 3]);


% 1. layer image
msk_mask = l3mask * 0.6;
img_mask = max(1 - l3mask, 0.4);
layimg = img_mask .* i1g + msk_mask .* sc(layermap, 'jet', [0, 3]);

% 2. obj image
objmap = A.object_map;
% TODO: assign random colour 
omask = objmap > 0;
o3mask = repmat(omask, [1, 1, 3]);
msk_mask = o3mask * 0.6;
img_mask = max(1 - o3mask, 0.4);
objimg = img_mask .* i1g + msk_mask .* sc(objmap, 'jet', [0, 3]);


%%% % %------------------------------------------------
%%% % % vasiliy color settings
%%% % %------------------------------------------------
%%% % jetmap = jet(4);
%%% % cyan = jetmap(2,:);
%%% % red = jetmap(3,:);
%%% % orange = jetmap(4,:);
%%% % green = [0.2,1.5,0.2];
%%% % red = [1 0 0]*1.2;
%%% % yellow = [1, 1, 0.0];
%%% % black = [0,0,0];
%%% % % cyan = red;
%%% % colormap{1} = [black; red; orange; yellow; cyan;  green];   % ours
%%% % colours = [colormap{1}; repmat(colormap{1},[50, 1])];
%%% % maxval = max(objmap(:));
%%% % seg = sc(double(objmap), colours(1:maxval+1,:) );
%%% % % I = repmat(im2double(rgb2gray(I)),[1 1 3]);
%%% % objimg2 = img_mask .* i1g + msk_mask .* seg;
%%% % keyboard;
%%% 
%%% % objim = cat(3, objmap == 1, objmap == 5, objmap == 16);
%%% % objim = 0.5 * im2single(objim) + o3mask * 0.5;
%%% % objim = sc(objim);
%%% % objimg2 = img_mask .* i1g + msk_mask .* objim;
%%% % keyboard;

% 3. fg prior img
MAXPFG = 2.0;
pfg = A.vis_prob_fg;
pfg_mask = pfg > 0;
pfg3 = repmat(pfg, [1, 1, 3]) * 0.6;
pfg3_mask = repmat(pfg_mask, [1, 1, 3]) * 0.6;
pfg_max = max(vec(pfg));
% pfgimg = max(1.0 - pfg3_mask, 0.4) .* i1g ...
%   + pfg3_mask .* sc(pfg, 'jet', [0, MAXPFG]);
pfgimg = max(1.0 - pfg3, 0.4) .* i1g ...
  + pfg3 .* sc(pfg, 'jet', [0, MAXPFG]);
pfgimg2 = sc(pfg, 'jet', [0, MAXPFG]);

% 4. bg prior img
MAXPBG = 2.0;
BG = isfield(A, 'vis_prob_bg');
if BG;
  pbg = A.vis_prob_bg;
  pbg_mask = pbg > 0;
  pbg3_mask = repmat(pbg_mask, [1, 1, 3]) * 0.6;
  pbg_max = max(vec(pbg));
  pbgimg = max(1.0 - pbg3_mask, 0.4) .* i1g ...
    + pbg3_mask .* sc(pbg, 'jet', [0, MAXPBG]);
end

% 5. cues img
MAXCON = 1.0; MINCON = 0.25; % MINCON = 0.025;
CON = isfield(A, 'constraints');
if CON;
  c = A.constraints;
  cw = A.constraint_weights;
  keep = cw > MINCON;
  c = c(keep, :);
  cw = cw(keep);
  
  cw = min(MAXCON, cw);
  % cw = ones(size(cw));
  [con1img, occ1img, occr1img] = vis_cues2(I1, c, cw, MINCON);
  i1 = sc(I1, 'gray');
  [con2img, occ2img, occr2img] = vis_cues2(i1, c, cw, MINCON);
end
if ~CON && FXF;
  conimg1 = A.vis_cues_btay;
end

% 6. weights img
wxy = vis_weights_occlusions(A.Wx, [], imsize, []);


% write our figure images
writepath = sprintf('./out/%s', seq);
% writepath = '/plot/btay/projects/detachable/cdov/figs';
outname = sprintf('%s_%s_%%s.png', ustr, fname);
if ~exist(writepath, 'dir'); mkdir(writepath); end;


o_img_file = fullfile(writepath, sprintf(outname, 'img'));
imwrite(i1g, o_img_file, 'png');
o_lay_file = fullfile(writepath, sprintf(outname, 'lay'));
imwrite(layimg , o_lay_file, 'png');
o_obj_file = fullfile(writepath, sprintf(outname, 'obj'));
imwrite(objimg , o_obj_file, 'png');
o_pfg_file = fullfile(writepath, sprintf(outname, 'pfg'));
imwrite(pfgimg, o_pfg_file, 'png');
o_pfg2_file = fullfile(writepath, sprintf(outname, 'pfg2'));
imwrite(pfgimg2, o_pfg2_file, 'png');
o_wxy_file = fullfile(writepath, sprintf(outname, 'wxy'));
imwrite(wxy, o_wxy_file, 'png');
if BG;
  o_pbg_file = fullfile(writepath, sprintf(outname, 'pbg'));
  imwrite(pbgimg, o_pbg_file, 'png');
end
if CON || FXF;
  o_con1_file = fullfile(writepath, sprintf(outname, 'con'));
  o_occ1_file = fullfile(writepath, sprintf(outname, 'occ'));
  o_occr1_file = fullfile(writepath, sprintf(outname, 'occr'));
  imwrite(con1img, o_con1_file, 'png');
  imwrite(occ1img, o_occ1_file, 'png');
  imwrite(occr1img, o_occr1_file, 'png');
  if CON;
    o_con2_file = fullfile(writepath, sprintf(outname, 'cong'));
    o_occ2_file = fullfile(writepath, sprintf(outname, 'occg'));
    o_occr2_file = fullfile(writepath, sprintf(outname, 'occrg'));
    imwrite(con2img, o_con2_file, 'png');
    imwrite(occ2img, o_occ2_file, 'png');
    imwrite(occr2img, o_occr2_file, 'png');
  end
end
end
