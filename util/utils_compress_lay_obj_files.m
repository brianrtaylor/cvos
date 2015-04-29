% load all files corresponding to this sequence
function lay = utils_compress_lay_files(dopath, seq, flag, rev_flag, index)
p = [dopath '/' seq '/'];
files = [];
if isempty(files)
  % try another naming style:
  files = dir( [p '/' seq '_*lay*_cdov_*', flag, '*.mat'] );
end
if isempty(files)
  files = dir( [p '/' seq '_*lay_*dov_frame*', flag, '*.mat'] );
end
lay = [];
T = length(files);
T2 = get_sequence_length(seq);
if (T < T2)
  fprintf('on %s : expected to see %d files, but only saw: %d\n', ...
    seq, T2, T);
  fprintf('pointless to continue.\n');
  return;
end
if nargin>=4 && strcmpi(rev_flag,'reverse')
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

if ~exist('index','var')
  for ii=1:T2
    fprintf('%s\n', files(ii).name);
    env = load([p, files(ii).name],'lay');
    if ii==1, lay = zeros( size(env.lay,1), size(env.lay,2), T2 ,'uint8'); end;
    lay(:,:,ii) = single(env.lay);
  end
else
  for ii=1:length(index)
    fprintf('%s\n', files( index(ii) ).name);
    env = load([p, files( index(ii) ).name],'lay');
    if ii==1, lay = zeros( size(env.lay,1), size(env.lay,2), length(index) ,'uint16'); end;
    lay(:,:,ii) = single(env.lay);
  end
end
end
