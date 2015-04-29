%-------------------------------------------------------------------------
% goal: to warp frame 2 output to frame 1 and frame N output to
% frame N - 1, thus we have segmentations on the whole sequence
%
% Use configuration for 2nd frame...
% @param: layers: frame 2 layers
% @param: uv: frame 1 uvf, warps 1 to 2
% @param: occ: frame 1 occf, warps 1 to 2
%
% ...or configuration for N-1th frame:
% @param: layers: frame N-1 layers
% @param: uv: frame N uvb, warps N to N-1
% @param: occ: frame N occb, warps N to N-1
%-------------------------------------------------------------------------
function layers_t0 = warp_output_layers(layers, uv, occ)
layers_t0 = utils_warp_image(layers, uv);
layers_t0 = remove_occb_mask(layers_t0, occ, 0, 0);
end