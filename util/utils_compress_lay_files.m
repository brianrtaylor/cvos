%-----------------------------------------------------------------------------
% utils_compress_lay_files
%
% aggregates layer and object segmentation results in individual .mat files 
% into a single file for evaluation 
%
% @return: lay (MxNxZ): stack of layer segmentations for Z sequence images
% @return: obj (MxNxZ): stack of object segmentations for Z sequence images
%
% @param: p: path to output files
% @param: seq: sequence name or prefix to output files
% @param: post: what string to append to the video name
% @param: rev_flag: whether or not to make a forward-backward noncausal video
% @param: index: which images to load (defaults: whole sequence)
%-----------------------------------------------------------------------------
function [lay, obj] = utils_compress_lay_files(p, seq, post, rev_flag, index)
str_search = fullfile(p, sprintf('%s_*_lay_cvos_%s.mat', seq, post));
files = dir(str_search);
lay = [];
T = length(files);
T2 = get_sequence_length(seq);
if (T < T2)
  fprintf('on %s : expected %d files, but only saw: %d\n', ...
    seq, T2, T);
  fprintf('pointless to continue.\n');
  keyboard;
  return;
end
if nargin>=4 && rev_flag
  if T ~= 2*T2-1
    error('not enoguh files.\n');
  end
  files = files(end:-1:1);
  files = files(1:T2);
end
if (T2 ~= T)
  fprintf('%s: expected to see %d files, but instead see %d\n', seq, T2, T);
  T = T2;
end

if ~exist('index','var');
  for ii=1:T2;
    fprintf('%s\n', files(ii).name);
    env = load(fullfile(p, files(ii).name), 'lay', 'object_map');
    if ii==1;
      lay = zeros(size(env.lay,1), size(env.lay,2), T2 ,'uint8'); 
      obj = zeros(size(env.object_map,1), size(env.object_map,2), ...
        T2 ,'uint16');
    end;
    lay(:,:,ii) = uint8(env.lay);
    obj(:,:,ii) = uint16(env.object_map);
  end
else
  for ii=1:length(index);
    fprintf('%s\n', files( index(ii) ).name);
    env = load(fullfile(p, files( index(ii) ).name), 'lay', 'object_map');
    if ii==1;
      lay = zeros(size(env.lay,1), size(env.lay,2), length(index) ,'uint8');
      obj = zeros(size(env.object_map,1), size(env.object_map,2), ...
        length(index) ,'uint16');
    end
    lay(:,:,ii) = uint8(env.lay);
    obj(:,:,ii) = uint16(env.object_map);
  end
end

% post process object map
obj1 = obj(:,:,1); obj1_keep = obj1;
vobj1 = unique(obj1(:,:,1));
vobj2 = unique(obj(:,:,2));
last_obj_id = max(vec(obj(:,:,end)));
noids = setdiff(vobj1, vobj2);
for k = 1:length(noids);
  last_obj_id = last_obj_id + 1;
  obj1(obj1_keep == noids(k)) = last_obj_id;
end
obj(:,:,1) = obj1;

if rev_flag;
  lay_file = fullfile(p, sprintf('%s_fb_lay.mat', seq));
  obj_file = fullfile(p, sprintf('%s_fb_obj.mat', seq));
else
  lay_file = fullfile(p, sprintf('%s_lay.mat', seq));
  obj_file = fullfile(p, sprintf('%s_obj.mat', seq));
end
save(lay_file, 'lay'); fprintf('%s file saved\n', lay_file);
save(obj_file, 'obj'); fprintf('%s file saved\n', obj_file);
end
