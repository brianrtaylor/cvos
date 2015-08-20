% function to draw avrious diagnostics of the boxes to see what's
% going on inside of them
%
% @note: assumes all box insides are the same size for now: 
%   prob, mask, etc.
function o = diagboxes(boxes)

% parse input
n_boxes = size(boxes, 1);
a = boxes(1).fg_mask;
[height, width] = size(a);
r = (height - 1) / 2;

% dividers
vbar1 = zeros(height, 1, n_boxes); vbar2 = vbar1;
vbar1(1:2:end) = 1;
vbar2(2:2:end) = 1;
vbar = [vbar1, vbar2];

hbar1 = zeros(1, 2 * width + 2, n_boxes); hbar2 = hbar1;
hbar1(1:2:end) = 1;
hbar2(2:2:end) = 1;
hbar = [hbar1; hbar2];

% output images
fg_mask = cat(3, boxes.fg_mask);
bg_mask = cat(3, boxes.bg_mask);
n_masks = size(fg_mask, 3);

fg_prob = cat(3, boxes.fg_prob);
bg_prob = cat(3, boxes.bg_prob);
n_probs = size(fg_prob, 3);

% colour
conf_colour = repmat(cat(3, boxes.conf_colour), height, width);
fg_prob_colour = cat(3, boxes.fg_prob_colour) .* conf_colour;
bg_prob_colour = cat(3, boxes.bg_prob_colour) .* conf_colour;
n_probs_colour = size(fg_prob_colour, 3);
if n_probs_colour < n_probs;
  pad0 = zeros(height, width, n_probs - n_probs_colour);
  fg_prob_colour = cat(3, fg_prob_colour, pad0);
  bg_prob_colour = cat(3, bg_prob_colour, pad0);
end

% shape
conf_shape = cat(3, boxes.conf_shape);
fg_prob_shape = cat(3, boxes.fg_prob_shape) .* conf_shape;
bg_prob_shape = cat(3, boxes.bg_prob_shape) .* conf_shape;
n_probs_shape = size(fg_prob_shape, 3);
if n_probs_shape < n_probs;
  pad0 = zeros(height, width, n_probs - n_probs_shape);
  fg_prob_shape = cat(3, fg_prob_shape, pad0);
  bg_prob_shape = cat(3, bg_prob_shape, pad0);
end

% put together output
omask = [fg_mask, vbar, bg_mask];
oprob = [fg_prob, vbar, bg_prob];
oprob_colour = [fg_prob_colour, vbar, bg_prob_colour];
oprob_shape = [fg_prob_shape, vbar, bg_prob_shape];

% old mask
% prob
% prob_colour
% prob_shape
% o = [omask; hbar; oprob_colour; hbar; oprob_uv; hbar; oprob_shape; hbar; oprob];
o = [omask; hbar; oprob_colour; hbar; oprob_shape; hbar; oprob];
end
