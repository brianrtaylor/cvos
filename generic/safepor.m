% performs probabilistic or on a and b
% C = A + B - A * B
function C = safepor(A,B)
% clean it
As = A;
Bs = B;
As(isnan(A)) = 0; % basically doesn't incorporate it
Bs(isnan(B)) = 0; 
C = (As + Bs) - (As .* Bs);
% C(isnan(A) & isnan(B)) = NaN; %? should I do this?
end
