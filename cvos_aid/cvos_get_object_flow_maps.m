%------------------------------------------------------------------
% objects top pre-processing for mean flow images
%------------------------------------------------------------------
function uvf_map = cvos_get_object_flow_maps(objects, object_map)
[rows, cols, ~] = size(object_map);
imsize = [rows, cols];

nobjects = max(object_map(:));
object_mean_uf    = nan(imsize);
object_mean_vf    = nan(imsize);
if nobjects > 0;
  object_ids = cat(1, objects.id);

  % counts by object ids and not just the order in the list
  for ok = 1:nobjects;
    obj = objects(object_ids == ok);
    object_mean_uf(object_map == ok) = obj.u;
    object_mean_vf(object_map == ok) = obj.v;
  end
end
uvf_map = cat(3, object_mean_uf, object_mean_vf);
end