function [files, valid] = purge_files(files, str)
N = length(files);
valid = true(N, 1);
for k = 1:N;
  if ~isempty(strfind(files(k).name, str));
    valid(k) = false;
  end
end
files = files(valid);
end