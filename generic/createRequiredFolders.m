function createRequiredFolders(d, decOpt)

% function will recursively required all necessary directories and sub-
% directories to build the folder d. d should be a specific filename
 
% turns periods to "p"s if set to true
if ~exist('decOpt', 'var');
  decOpt = false;
end

[p, name, ext] = fileparts(d);
if ~exist(p, 'dir') && ~isempty(p);
  p
  createRequiredFolders(p);
end

if isempty(ext);
  if decOpt;
    d = strrep(d, '.', 'p');
  end
  mkdir(d);
  fprintf('%s created\n', d);
end
