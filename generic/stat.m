function stat(A);
% just reports some stats, most useful is max, min

out = struct;
Avec = double(A(:));

out.max = max(Avec);
out.min = min(Avec);
out.mean = mean(Avec);
out.variance = var(Avec);

out.elements = unique(Avec);
out.numUniqueElements = length(out.elements);
out.dims = size(A);
out.numElements = length(Avec);

ndims = length(out.dims);

if ndims == 2;
  fprintf('dims: %d x %d, ', out.dims);
else ndims == 3;
  fprintf('dims: %d x %d x %d, ', out.dims);
end

fprintf('(%d total elements)\n', out.numElements);
fprintf('min: %8.6f,  max: %8.6f\n', out.min, out.max);
fprintf('mean: %f, variance: %f\n', out.mean, out.variance);
fprintf('%d unique elements: ', out.numUniqueElements);

numToShow = min(out.numUniqueElements, 10);
fprintf('%d, ', out.elements(1:numToShow - 1));
fprintf('%8.4f\n', out.elements(numToShow));
