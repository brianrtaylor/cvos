%---------------------------------------------------------------------
% expands occlusion region to region likely needing flow fixing
%
% @note: solely for cross bilateral filter step
%---------------------------------------------------------------------
function region = occ_expand(occ, uv)
[M, N] = size(occ);
mag_nuv = sqrt(sum(uv .^ 2, 3));
region = occ;
for y = 1:M;
  for x = 1:N;
    r = occ(y, x);
    Q = max(1, floor(mag_nuv(y, x) * r));
    Q = clip(Q, 2, 50);
    ymin = max(1, y-Q);
    xmin = max(1, x-Q);
    ymax = min(M, y+Q);
    xmax = min(N, x+Q);
    region(y, x) = max(vec(occ(ymin:ymax, xmin:xmax)));
  end
end
end
