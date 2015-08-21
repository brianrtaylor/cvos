% convenient struct definition for storing past information
function past = struct_past
past = struct;
past.uvf_cbf = [];
past.uvf = [];
past.uvf_rev = [];

past.uvb = [];
past.uvb_rev = [];
past.uvf = [];
past.uvf_rev = [];

past.occf_cbf = [];
past.occf_rev = [];

past.constraints_b = [];
past.constraint_weights_b = [];
past.constraints_causal_b = [];
past.constraint_weights_causal_b = [];
past.constraints_f = [];
past.constraint_weights_f = [];
past.constraints_causal_f = [];
past.constraint_weights_causal_f = [];

past.weights = [];
past.constraints = [];
past.constraint_weights = [];
past.constraint_ages = [];
past.prob_fg = [];
past.xi_b = [];
past.xi_f = [];
past.layers = [];
past.prob_fg = [];
past.prob_fg_layer = [];

past.unity = [];
past.unity_layer = [];

% things warped to current domain in past, are called past.t0.
past.t0.layers = [];
past.t0.weights = [];
end
