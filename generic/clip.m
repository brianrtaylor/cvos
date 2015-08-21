%-----------------------------------------------------------------------------
% clip
%
% clips values in the input to a specific range
%
% @return: B: output clamped matrix
% @param: A: input matrix
% @param: lo: low clamp value
% @param: hi: high clamp value
%-----------------------------------------------------------------------------
function B = clip(A, lo, hi)
if ~exist('lo', 'var'); lo = 0.0; end;
if ~exist('hi', 'var'); hi = 1.0; end;
B = min(hi, max(A, lo));
end
