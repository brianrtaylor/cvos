%	computes the sigmoid if the input x given offset and scale parameters
function y = sigmoid(x, offset, scale)
y = 1 - 1 ./(1 + exp(-scale .* (offset - x)));
end
