% focuses on removing not as much for smaller regions, whereas large 
% occb regions have more influence
%
% @goal: to prevent small occlusion cues in the middle of an object 
%   to hurt foreground prior and same layer weight
function [out, out_msfm] = remove_occb_mask(in, occb, erodesz, NOCEIL)

if ~exist('NOCEIL', 'var'); NOCEIL = false; end;

if ~exist('erodesz', 'var'); erodesz = 0; end;
[rows, cols] = size(occb);
npx = rows * cols;
occb_mask = occb > 1e-6;

b = in;
b(isnan(b)) = 0;
bb = b > 1e-6;
[ii, jj] = find(1 - bb);
bb = msfm(double(bb), [ii, jj]', true, true);

[ii, jj] = find(1 - occb_mask);
ob = msfm(double(occb_mask), [ii, jj]', true, true);

obcc = bwlabel(occb_mask);

% tic
% ob1 = zeros(size(ob));
% for ok = 1:max(obcc(:));
%   maxx = max(vec(ob(obcc == ok)));
%   ob1(obcc == ok) = maxx;
% end
% 
% bbb = bb;
% % bbb(occb_mask) = max(0.0, bbb(occb_mask) - (1 + 2 * ob1(occb_mask)));
% bbb = max(0.0, bbb - (1 + 2 * ob1));
% toc

% tic
[bbb, occ_cc_max_mex] = remove_occb_mask_loop_mex( ...
  double(bb), double(ob), double(obcc));
% toc

% % [out_msfm, obcc_maxe] = fxn(occ_msfm (ob), in_msfm (bbb))
% tic
% obcc_maxe = zeros(max(obcc(:)) + 1, 1);
% for ii = 1:npx;
%   l = obcc(ii) + 1;
%   obcc_maxe(l) = max(obcc_maxe(l), ob(ii));
% end
% 
% bbb = bb;
% % for ii = 1:npx;
% %   bbb(ii) = bbb(ii) - (1 + 2 * obcc_maxe(obcc(ii) + 1));
% % end
% ind = (obcc > 0);
% bbb(ind) = bb(ind) - (1 + 2 * obcc_maxe(obcc(ind) + 1));
% toc

% bbb = medfilt2(ceil((bbb > 0.0) .* b - 0.1), [3, 3]);
if NOCEIL;
  out = (bbb > erodesz) .* b;
else
  out = ceil((bbb > erodesz) .* b - 0.1);
end
out_msfm = bbb;
end