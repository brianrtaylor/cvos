% creates a figure without stealing focus from the matlab command window
function h = fig(varargin)
if mod(numel(varargin),2)
  f = varargin{1};
  parameters = varargin(2:end);
else
  parameters = varargin;
  f = [];
end

if ishandle(f)
  set(0,'CurrentFigure',f);
elseif ~isempty(f)
  f = figure(f);
else
  f = figure;
end

if ~isempty(parameters)
  set(f,parameters{:});
end

if nargout == 1
  h = f;
end
end
