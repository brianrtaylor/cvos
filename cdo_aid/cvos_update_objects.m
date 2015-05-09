%-------------------------------------------------------------------------
% cvos_update_objects
%
% % @param: objects : structure array of objects
% @param: object_map : image of object labels obtained from 
%   layers_to_detachable_objects fxn call
% @param: uvf : flow from t to t+1
%
% @note: object representation
% * id : id of object
% * u, v : average motion of the object
%
%
% @problem: what if we have objects, they disappear later, and then lastly
% we get objects again. what will happen? --> they become different objects
%-------------------------------------------------------------------------
function [objects_out, object_map_out, n_active] = cvos_update_objects( ...
  objects, object_map, object_map_snap, uvb_warp, uvf)
nnewobjects = max(object_map_snap(:));
% nobjects = length(objects);

ids_alive_objects = vec(unique(object_map(:)));
if ids_alive_objects(1) == 0; ids_alive_objects(1) = []; end;
if isempty(ids_alive_objects); ids_alive_objects = []; end;
n_alive_objects = length(ids_alive_objects);

[rows, cols, ~] = size(uvf);
imsize = [rows, cols];

u = uvf(:, :, 1);
v = uvf(:, :, 2);

%-------------------------------------------------------------------------
% warp old object map from past frame into frame t0
%-------------------------------------------------------------------------
object_map_t0 = utils_warp_image(double(object_map), uvb_warp);
object_map_t0(isnan(object_map_t0) | isinf(object_map_t0)) = 0.0;
object_map_out = zeros(imsize, 'single');

%-------------------------------------------------------------------------
% for each object in object_map_snap, attach object or create new one
%-------------------------------------------------------------------------
ids = []; max_id = 0;
if ~isempty(objects);
  ids = cat(1, objects.id); max_id = max(ids);
end

% lifted representation of alive object maps
obj_maps = false([rows, cols, n_alive_objects]);
for j = 1:n_alive_objects;
  obj_maps(:, :, j) = (object_map_t0 == ids_alive_objects(j));
end

ids_new = []; nadd = 0;
ids_still_alive_objects = []; n_still_alive_objects = 0;
for k = 1:nnewobjects;
  % 1. assign an object to this snap_object or create a new object for it
  snap_mask = object_map_snap == k;
  npx = sum(snap_mask(:));
  
  % find object with greatest overlap
  f1max = -1.0;

  if n_alive_objects > 0;
    f1 = zeros(n_alive_objects, 1, 'single');
    for j = 1:n_alive_objects;
      dat = evaOne(snap_mask, obj_maps(:, :, j));
      f1(j) = dat.f;
    end
    f1(isnan(f1)) = -1.0;
    [f1max, id_ind] = max(f1);
  end

  
  % decide if old object (high overlap) or new one (low overlap)
  if (f1max > 0.5); % old object
    snap_id = ids_alive_objects(id_ind);

    ids_still_alive_objects = cat(1, ids_still_alive_objects, snap_id);
    n_still_alive_objects = n_still_alive_objects + 1;
  else % new object
    snap_id = max_id + 1;
    max_id = snap_id;

    ids_new = cat(1, ids_new, snap_id);
    nadd = nadd + 1;
    % ids = cat(1, ids, snap_id);
    % ids_alive_objects_t0 = cat(1, ids_alive_objects_t0, snap_id);
  end

  % fill in object into output object map
  object_map_out(snap_mask) = snap_id;
end
% nobjects = nobjects + nadd;
ids_alive_objects_t0 = cat(1, ids_alive_objects, ids_new);
n_alive_objects_t0 = n_alive_objects + nadd;

ids_spawn_objects = ids_new;
n_spawn_objects = nadd;

%-------------------------------------------------------------------------
% post-process output object map to make sure no disconnected componenets
% have the same object id
%-------------------------------------------------------------------------
object_map_to_add = zeros(imsize, 'single');
nadd = 0; ids_new = [];
for j = 1:n_alive_objects_t0;
  [L, num] = bwlabel(object_map_out == ids_alive_objects_t0(j));
  if num > 1; % split to create a new object
    for jj = 1:num;
      snap_id = max_id + 1;
      max_id = snap_id;
      nadd = nadd + 1;
      % ids = cat(1, ids, snap_id);
      ids_new = cat(1, ids_new, snap_id);
      object_map_to_add(L == jj) = snap_id;
    end
  end
end
ind = object_map_to_add > 0;
object_map_out(ind) = object_map_to_add(ind);

ids_spawn_objects = cat(1, ids_spawn_objects, ids_new);
n_spawn_objects = n_spawn_objects + nadd;

%-------------------------------------------------------------------------
% update objects after fully deciding on new object map
%-------------------------------------------------------------------------
% update flow for persisting objects
if n_still_alive_objects > 0;
  ids = cat(1, objects.id);
  % 1. get indices
  % TODO: can make this faster
  for j = 1:n_still_alive_objects;
    ind_objects = find(ids == ids_still_alive_objects(j));

    mask = object_map_out == ids_still_alive_objects(j);
    npx = sum(mask(:));

    ou = u(mask);
    ov = v(mask);
    
    mou = sum(ou(:)) / npx;
    mov = sum(ov(:)) / npx;

    objects(ind_objects).u = single(mou);
    objects(ind_objects).v = single(mov);
  end
end

% 2. create and compute flow for new objects
objects_spawn = [];
for j = 1:n_spawn_objects;
  mask = object_map_out == ids_spawn_objects(j);
  npx = sum(mask(:));
  
  ou = u(mask);
  ov = v(mask);
  
  mou = sum(ou(:)) / npx;
  mov = sum(ov(:)) / npx;
 
  obj_new = struct;
  obj_new.id = single(ids_spawn_objects(j));
  obj_new.u = single(mou);
  obj_new.v = single(mov);
  
  objects_spawn = cat(1, objects_spawn, obj_new);
end

objects_out = cat(1, objects, objects_spawn);
n_active = n_still_alive_objects + n_spawn_objects;
end
