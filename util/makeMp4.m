%-----------------------------------------------------------------------------
% small helper function to make your movies
% example: avconv -r 5 -f image2 -i cars5_%03d_causal_0.png cars5_causal_0.mp4
%-----------------------------------------------------------------------------
function makeMp4(imgpth, seq, endname, framerate, bitrate)

if ~exist('framerate', 'var'); framerate= 5; end;
if ~exist('bitrate', 'var'); bitrate = framerate * 500; end;

search_str = fullfile(imgpth, sprintf('%s_*_%s.png', seq, endname));
fs = dir(search_str);
N = length(fs);
ndigits = ceil(log10(N)) + 1;

cmd = '';
tmp_str = sprintf('.%%0%dd.png', ndigits);
mp4_file_name = sprintf('%s_%s_%02dfps.mp4', seq, endname, framerate);
cmd = sprintf('avconv -r %d -f image2 -i %s -b %dk %s', ...
  framerate, tmp_str, bitrate, mp4_file_name);
% cmd = sprintf('avconv -r %d -f image2 -i %s %s', ...
%   framerate, tmp_str, mp4_file_name); % works

% setup commands to rename the files to .00{1,2,3,...}.png and back
mvcmd = '';
mvbackcmd = '';
for k = 1:N;
  name1 = fullfile(imgpth, fs(k).name);
  name2 = fullfile(imgpth, sprintf(tmp_str, k));
  mvcmd     = [mvcmd,     sprintf('mv %s %s; ', name1, name2)];
  mvbackcmd = [mvbackcmd, sprintf('mv %s %s; ', name2, name1)];
end

% rename
unix(mvcmd);

here = pwd();
cd(imgpth);
if exist(mp4_file_name, 'file'); delete(fullfile(imgpth, mp4_file_name)); end;
atic = tic();
unix(cmd);
atic = toc(atic);
cd(here);

% rename back
unix(mvbackcmd);
fprintf('%s created in %0.3f s. Enjoy!\n', mp4_file_name, atic);
end

