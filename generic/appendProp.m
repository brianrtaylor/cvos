function [s, ssub] = appendProp(strin, val, numtype, nsigfigs)
s = ''; ssub = '';
if ~exist('nsigfigs', 'var'); nsigfigs = 3; end;
if ~exist('type', 'var'); numtype = 'float'; end;

if isstr(val);
  ssub = val;
elseif isnumeric(val); %TODO: negative values
  if strcmp(numtype, 'int') == 1; % not float
    if (floor(val) ~= val); 
      fprintf('appendProp: val not a positive integer input. val = %f\n', val);
      keyboard;
      return;
    end
    ssub = sprintf(sprintf('%%0%dd', nsigfigs), val);
  else
    ssub = strrep(sprintf(sprintf('%%0.%df', nsigfigs), val), '.', 'p');
  end
else;
  fprintf('sorry, meaningless input');
  keyboard;
end
s = sprintf('%s%s', strin, ssub);
end
