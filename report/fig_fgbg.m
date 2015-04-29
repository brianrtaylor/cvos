function fig_gfbg(seq, id, ustr, outpath)
if ~exist('ustr', 'var'); ustr = 'check'; end;
if ~exist('outpath', 'var'); outpath = '/plot/btay/projects/detachable/cdov/week20-cdov-r0.6-fullbg-moseg/'; end;

nameStr = '%s_%06d';
fname = sprintf(nameStr, seq, id);

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

f = dir(fullfile(outpath, seq, sprintf('%s_*_lay_*.mat', seq)));
mat_file = fullfile(outpath, seq, f(id).name);
A = load(mat_file);

%-------------------------------------------------------------------------
% make images 
%-------------------------------------------------------------------------
% layer image
layermap = A.lay;
layermap(1,1,:)=0;
layermap(1,2,:)=3;
layermap(1:10,1:10)=0;
layermap(1:10,11:20)=1;
layermap(1:10,21:30)=2;
layermap(1:10,31:40)=3;
layermap = min(3, layermap); % cuts off anything higher
lmask = (layermap > 0);
l3mask = repmat(lmask, [1, 1, 3]);


% 1. layer image
msk_mask = l3mask * 0.6;
img_mask = max(1 - l3mask, 0.4);
layimg = img_mask .* i1g + msk_mask .* sc(layermap, 'jet');

% 2. obj image
objmap = A.object_map;
% TODO: assign random colour 
omask = objmap > 0;
o3mask = repmat(omask, [1, 1, 3]);
msk_mask = o3mask * 0.6;
img_mask = max(1 - o3mask, 0.4);
objimg = img_mask .* i1g + msk_mask .* sc(objmap, 'jet');

% 3. fg prior img
MAXPFG = 1.0;
pfg = A.vis_prob_fg;
pfg_mask = pfg > 0;
pfg3 = repmat(pfg, [1, 1, 3]) * 0.6;
pfg3_mask = repmat(pfg_mask, [1, 1, 3]) * 0.6;
pfg_max = max(vec(pfg));
% pfgimg = max(1.0 - pfg3_mask, 0.4) .* i1g ...
%   + pfg3_mask .* sc(pfg, 'jet', [0, MAXPFG]);
pfgimg = max(1.0 - pfg3, 0.4) .* i1g ...
  + pfg3 .* sc(pfg, 'jet', [0, MAXPFG]);

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



% write our figure images
writepath = './out';
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
if BG;
  o_pbg_file = fullfile(writepath, sprintf(outname, 'pbg'));
  imwrite(pbgimg, o_pbg_file, 'png');
end
end
