% @param: layers: integer map larger than 0

function [object_map, last_obj_id] = layers_to_detachable_objects(layers)
if nargin == 0;
  test_layers_to_detachable_objects; return;
end

[rows, cols, ~] = size(layers);
object_map = zeros(rows, cols);
max_layer = max(layers(:));
last_obj_id = 0;
for k = 1:max_layer;
  mask_k = abs(layers - k) < 0.5;
  [lay_cc, num_cc] = bwlabel(mask_k);
  object_map = object_map + (lay_cc + mask_k * last_obj_id); 
  last_obj_id = last_obj_id + num_cc;
end
end


%-------------------------------------------------------------------------
% function test function
%-------------------------------------------------------------------------
function test_layers_to_detachable_objects

% make test image
A = zeros(100, 100);
b1y = 20:70; b1x = 10:40;
b2y = 40:60; b2x = 30:60;
b3y = 70:90; b3x = 70:90;
A(b1y, b1x) = 1;
A(b2y, b2x) = 2;
A(b3y, b3x) = 1;

% try it
object_map = layers_to_detachable_objects(A);

% check
assert(all(vec(object_map(b2y, b2x) == 3)), 'box 2 messed up');
assert(all(vec(object_map(b3y, b3x) == 2)), 'box 3 messed up');
assert(all(vec(object_map(b1y, b1x) >= 1)), 'box 1 messed up');

% win
fprintf('%s: all tests passed\n', mfilename);
end
