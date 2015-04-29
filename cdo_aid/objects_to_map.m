%
% @param: objects:
% * id : id in sequence
% * mask : where object was (binary mask)
% * u : mean u
% * v : mean v

function [object_map] = objects_to_map(objects, imsize)
nobjects = size(objects, 1);
object_map = zeros(imsize);

[~, lay_ind] = sort(cat(1, objects.lay), 1, 'ascend');
objects = objects(lay_ind, 1);

for k = 1:nobjects;
  o = objects(k);

  object_map(o.mask) = o.id;
end
end