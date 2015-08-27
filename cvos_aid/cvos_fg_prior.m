%-----------------------------------------------------------------------------
% cvos_fg_prior
%
% computes the prior on which pixels in the image domain are foreground based
% on the foreground region from the prior frame
%
% @return: prob_fg: the foreground probability map over the image domain
% @param: past: struct holding information from prior frames
% @param: uvb: optical flow from frame t to t-1
% @param: occb_mask: the region unoccluded when moving from frame t-1 to t
% @param: opts : options and parameters
%   * OCC_INSIDE_FG_MEMORY_FACTOR
%   * FOREGROUND_WARP_WEIGHT
%   * FG_MEMORY_FACTOR
%   * PROB_FOREGROUNDA
%   * CAUSAL
%   * dsk
%   * VIS
%-----------------------------------------------------------------------------
function prob_fg = cvos_fg_prior(past, uvb, occb_mask, opts)
v2struct(opts);

% parameters
ERODE_SZ = 3.0;
g_tiny = fspecial('gaussian', 5, 1);
mag_uvb_squared = sum(uvb .^ 2, 3);
mag_uvb = sqrt(mag_uvb_squared);
pfg = past.prob_fg;

%-----------------------------------------------------------------------------
% diminish prob_fg (age, warp)
%-----------------------------------------------------------------------------
if DO_FOREGROUND_WARP_WEIGHT;
  pfg = past.prob_fg .* past.w_fg_uvf .* past.w_warp_uvf;
end
pfg = pfg * (1 - FG_MEMORY_FACTOR);

%-----------------------------------------------------------------------------
% boost prob_fg if correct
%-----------------------------------------------------------------------------
pfg = pfg + past.prob_fg_layer;

%-----------------------------------------------------------------------------
% warp to current frame, the remove occb + erode
%-----------------------------------------------------------------------------
fg_warped = utils_warp_image(pfg, uvb);
fg_warped(isnan(fg_warped)) = 0;

% remove occb
unity_mask = clip(past.t0.layers, 0.0, 1.0);
occb_mask = max(0.0, occb_mask - unity_mask * OCC_INSIDE_FG_MEMORY_FACTOR);
[fg_warped_1, fg_msfm] = remove_occb_mask(fg_warped, occb_mask, ERODE_SZ, 1);

assert(all(vec(fg_warped  ) >= 0),   'fg_warped < 0');
assert(all(vec(fg_warped_1) >= 0), 'fg_warped_1 < 0');

if VIS < 150;
  fig(115); clf; imagesc([fg_warped, fg_warped_1]);
  title('fg warped | 1'); colorbar; drawnow;
end

% erode
prob_fg = fg_warped_1;
prob_fg( isnan(prob_fg ) ) = 0;

% erode edges based on flow magnitude in that region
prob_fg_chk = imerode(prob_fg > 0.01, dsk);
diminish_fg = 1 - exp(- max(0.0, fg_msfm - 2 * ERODE_SZ) ...
  ./ max(1.0, (mag_uvb / 2.0)));
diminish_fg = diminish_fg .* prob_fg_chk;
prob_fg = prob_fg .* diminish_fg;
prob_fg = imfilter(prob_fg, g_tiny);
end
