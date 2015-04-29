%----------------------------------------------------------------------------%
% ttitle(txt, color, x, y)
%
% helper function to place title text in a vl_tightsubplot image
%
% @param: txt : required
% @param: color : what color, defaults to blue
% @param: x : x position, defaults to around middle
% @param: y : y position, defaults to 40 pixels from top
%----------------------------------------------------------------------------%
function ttitle(txt, color, x, y)
s = get(gcf, 'Position');
if ~exist('x', 'var'); x = s(3); end;
if ~exist('y', 'var'); y = 40; end;
if ~exist('color', 'var'); color = 'blue'; end;
% just hope they might only dump numbers here
if strcmp(class(txt), 'char') ~= 1; txt = num2str(txt); end;
h = text((x/2 - (length(txt) * 6) / 2), y, txt);
set(h, 'Color', color);
end
