function fig_cbf(seq, id, model, ustr, outpath)
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

f_file = sprintf('%s_*%06d*_lay_*%s*.mat', seq, id, model);
f_srch = fullfile(outpath, seq, f_file);
f = dir(f_srch);
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

flo = A.flo;
v2struct(flo); 
% flo.uvf           = single(uvf);
% flo.uvf_cbf       = single(uvf_cbf);
% flo.uvf_bflt      = single(uvf_bflt);
% flo.uvf_cbf_bflt  = single(uvf_cbf_bflt);
% 
% flo.uvb           = single(uvb);
% flo.uvb_cbf       = single(uvb_cbf);
% flo.uvb_bflt      = single(uvb_bflt);
% flo.uvb_cbf_bflt  = single(uvb_cbf_bflt);
% 
% flo.occf          = single(occf);
% flo.occf_cbf      = single(occf_cbf);
% flo.occf_cbf_prob = single(occf_cbf_prob);
% flo.occb          = single(occb);
% flo.occb_cbf      = single(occb_cbf);
% flo.occb_cbf_prob = single(occb_cbf_prob);


uv = flowToColor([[uvb, uvb_bflt, uvf, uvf_bflt]; ...
  [uvb_cbf, uvb_cbf_bflt, uvf_cbf, uvf_cbf_bflt]]);
uv1 = flowToColor([[uvf, uvf_bflt]; [uvf_cbf, uvf_cbf_bflt]]);
uv2 = flowToColor([[uvb, uvb_bflt]; [uvb_cbf, uvb_cbf_bflt]]);

% fig(2); imagesc(uv);
% fig(3); imagesc(uv1);
% fig(4); imagesc(uv2);


% return

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

% 5. cues img
MAXCON = 1.0; MINCON = 0.025;
CON = isfield(A, 'constraints');
if CON;
  c = A.constraints;
  cw = A.constraint_weights;
  keep = cw > MINCON;
  c = c(keep, :);
  cw = cw(keep);
  
  % cw = min(MAXCON, A.constraint_weights);
  cw = ones(size(cw));
  con1img = vis_cues2(I1, c, cw, MINCON);
  i1 = sc(I1, 'gray');
  con2img = vis_cues2(i1, c, cw, MINCON);
end
if ~CON && FXF;
  conimg1 = A.vis_cues_btay;
end
  


%------------------------------------------------------
% write our figure images
%------------------------------------------------------
writepath = sprintf('./out/%s', seq);
% writepath = '/plot/btay/projects/detachable/cdov/figs';
outname = sprintf('%s_%s_%%s.png', ustr, fname);
if ~exist(writepath, 'dir'); mkdir(writepath); end;


% writing flows
[rows, cols, ~] = size(uvf);
uvfp = flowToColor([uvf, uvf_cbf]);
uvf1 = uvfp(1:rows, 1:cols, :);
uvf2 = uvfp(1:rows, (cols+1):(2*cols), :);
o_flo_file = fullfile(writepath, sprintf(outname, 'uvf1'));
imwrite(uvf1, o_flo_file, 'png');
o_flo_file = fullfile(writepath, sprintf(outname, 'uvf2'));
imwrite(uvf2, o_flo_file, 'png');

o_occb_file = fullfile(writepath, sprintf(outname, 'occb'));
o_occf_file = fullfile(writepath, sprintf(outname, 'occf'));
% o_occbcbf_file = fullfile(writepath, sprintf(outname, 'occbcbf'));
% o_occfcbf_file = fullfile(writepath, sprintf(outname, 'occfcbf'));

img_eb = sc(occb_cbf_prob, 'jet');
img_ef = sc(occf_cbf_prob, 'jet');

imwrite(uint8( 255.0 * img_eb), o_occb_file, 'png');
imwrite(uint8( 255.0 * img_fg), o_occf_file, 'png');
% imwrite(occbcbf, o_occbcbf_file, 'png');
% imwrite(occfcbf, o_occfcbf_file, 'png');


% need to get constraints for both cases ehhh... gg, have to get non bf
% flow constraints later if this flows look good. add constraints change
% later, start with just flow improvements




% rest of old stuff
o_img_file = fullfile(writepath, sprintf(outname, 'img'));
imwrite(i1g, o_img_file, 'png');
o_lay_file = fullfile(writepath, sprintf(outname, 'lay'));
imwrite(layimg , o_lay_file, 'png');
o_obj_file = fullfile(writepath, sprintf(outname, 'obj'));
imwrite(objimg , o_obj_file, 'png');

if CON || FXF;
  o_con1_file = fullfile(writepath, sprintf(outname, 'con'));
  imwrite(con1img, o_con1_file, 'png');
  if CON;
    o_con2_file = fullfile(writepath, sprintf(outname, 'cong'));
    imwrite(con2img, o_con2_file, 'png');
  end
end

o_flo_file
end
