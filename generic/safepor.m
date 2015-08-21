% "safely" computes "probabilistic or" of inputs A and B (C = A + B - A * B)
% essentially cleans the input of nans
function C = safepor(A,B)
As = A;
Bs = B;
As(isnan(A)) = 0; % remove nans
Bs(isnan(B)) = 0; 
C = (As + Bs) - (As .* Bs);
end
