%-----------------------------------------------------------------------------
% make_difference_operator
%
% USAGE: [Dx, Dy, dx_inds, dy_inds] = make_difference_operator(imsize)
%
% @return: D{x,y}: sparse difference operators for x and y directions
% @return: d{x,y}_inds: indices for edges in x and y directions
% @param: imsize: image dimensions in (height, width)
%-----------------------------------------------------------------------------
function [Dx, Dy, dx_inds, dy_inds] = make_difference_operator(imsize)
x = 1:imsize(2);
y = 1:imsize(1);

[xx,yy] = meshgrid(x,y);
xx_ = circshift(xx,[0, -1]);
yy_ = circshift(yy,[-1, 0]);

dx1 = [xx(:), yy(:)];
dx2 = [xx_(:), yy(:)];

dy1 = [xx(:), yy(:)];
dy2 = [xx(:), yy_(:)];
    
dx1 = sub2ind(imsize(1:2), dx1(:, 2), dx1(:, 1));
dx2 = sub2ind(imsize(1:2), dx2(:, 2), dx2(:, 1));
dy1 = sub2ind(imsize(1:2), dy1(:, 2), dy1(:, 1));
dy2 = sub2ind(imsize(1:2), dy2(:, 2), dy2(:, 1));

m = numel(dx1);
dx_ = [dx1; dx2];
wx_ = [-ones(m, 1); ones(m, 1)];
dy_ = [dy1; dy2];
wy_ = [-ones(m, 1); ones(m, 1)];
zz = [(1:numel(dx1))'; (1:numel(dx1))'];

onevec = ones(m,1);
nedges = numel(dx1);
nnodes = imsize(1) * imsize(2);

Dx = sparse(zz, dx_, wx_, nedges, nnodes, 2 * nedges);
Dy = sparse(zz, dy_, wy_, nedges, nnodes, 2 * nedges);

dx_inds = [dx1, dx2];
dy_inds = [dy1, dy2];
end
