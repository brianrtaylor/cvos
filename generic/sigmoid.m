function y = sigmoid(x, offset,scale)
%		FUNCTION FX = SIGMOID(X, OFFSET, SCALE)
y = 1-1./(1 + exp(-scale.*(offset-x)));
end


%%% function fx = sigmoid(x, amplitude, slope, phase)
%%% %		FUNCTION FX = SIGMOID(X, AMPLITUDE, SLOPE, PHASE)
%%% 
%%% fx = 1 ./ (1 + phase * exp(-slope * (x)));
%%% fx = amplitude * fx;
