%-------------------------------------------------------------------------
% cvos_object_map_to_objects
%
% % @param: objects : structure array of objects
% @param: object_map : image of object labels obtained from 
%   layers_to_detachable_objects fxn call
% @param: uvf : flow from t to t+1
%
% @note: object representation
% * id : id of object
% * u, v : average motion of the object
%-------------------------------------------------------------------------
function objects = cvos_object_map_to_objects(object_map, uvf)
nobjects = max(object_map(:));
objects = [];

u = uvf(:, :, 1);
v = uvf(:, :, 2);

for k = 1:nobjects;
  mask = object_map == k;
  npx = sum(mask(:));
  
  ou = u(mask);
  ov = v(mask);
  
  mou = sum(ou(:)) / npx;
  mov = sum(ov(:)) / npx;
 
  obj_new = struct;
  obj_new.id = k;
  obj_new.u = mou;
  obj_new.v = mov;
  
  objects = cat(1, objects, obj_new);
end
end