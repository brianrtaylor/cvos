% computes "probabilistic or" of probabilities a and b (C = A + B - A * B)
function C = por(A,B)
C = (A + B) - (A .* B);
end
