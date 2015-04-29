%------------------------------------------------------------------
% objects top pre-processing for mean flow images
%------------------------------------------------------------------
function uvf_map = tao_get_object_flow_maps(objects, object_map)
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

% % % % function [uvf_map, uvf_map_t0] = tao_get_object_flow_maps( ...
% % % %   objects, object_map, uvb)
% % % % [rows, cols, ~] = size(object_map);
% % % % imsize = [rows, cols];
% % % % 
% % % % object_map_t0 = round(utils_warp_image(object_map, uvb));
% % % % object_map_t0(isnan(object_map_t0) | isinf(object_map_t0)) = 0;
% % % % nobjects = max(object_map_t0(:));
% % % % % object_mean_uf_t0 = zeros(imsize);
% % % % % object_mean_vf_t0 = zeros(imsize);
% % % % % object_mean_uf    = zeros(imsize);
% % % % % object_mean_vf    = zeros(imsize);
% % % % object_mean_uf_t0 = nan(imsize);
% % % % object_mean_vf_t0 = nan(imsize);
% % % % object_mean_uf    = nan(imsize);
% % % % object_mean_vf    = nan(imsize);
% % % % if nobjects > 0;
% % % %   object_ids = cat(1, objects.id);
% % % %   
% % % %   % counts by object ids and not just the order in the list
% % % %   for ok = 1:nobjects;
% % % %     obj = objects(object_ids == ok);
% % % %     object_mean_uf_t0(object_map_t0 == ok) = obj.u;
% % % %     object_mean_vf_t0(object_map_t0 == ok) = obj.v;
% % % %     object_mean_uf(object_map == ok) = obj.u;
% % % %     object_mean_vf(object_map == ok) = obj.v;
% % % %   end
% % % % end
% % % % uvf_map_t0 = cat(3, object_mean_uf_t0, object_mean_vf_t0);
% % % % uvf_map = cat(3, object_mean_uf, object_mean_vf);
% % % % end
