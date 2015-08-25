%-----------------------------------------------------------------------------
% aggregate_pairs_fast
%
% given a list of pairs (or edges), this function will combine repeats into 
% a single entry and aggregate their associated weights via mex functions
%
% @return: unique_pairs: the remaining unique pairs
% @return: score_pairs: the weights associated with the remaining unique pairs
% @return: ind_pairs: indices corresponding to the elements of unique_pairs
% @param: pairs: edges between node 1 (1st column) and node 2 (2nd column)
% @param: pair_weights: weights associated with each edge
%-----------------------------------------------------------------------------
function [unique_pairs, score_pairs, ind_pairs] = aggregate_pairs_fast( ...
  pairs, pair_weights)

if isempty(pairs) && isempty(pair_weights)
  unique_pairs = pairs;
  score_pairs = pair_weights;
  ind_pairs = [];
  return
end

% test function 
if nargin == 0; test_aggregate_pairs; return; end;
npairs = size(pairs,1);

if ~exist('pair_weights', 'var');
    pair_weights = ones(npairs, 1);
end;
[pairs, sorti] = sortrows(pairs);
pair_weights = pair_weights(sorti);
unique_pairs = unique(pairs, 'rows');

nunique_pairs = size(unique_pairs,1);

[unique_pairs, score_pairs, ind_pairs] = aggregate_pairs_loop_mex( ...
  pairs, unique_pairs, pair_weights, sorti);
fprintf('%s: pruned edges: %d --> %d\n', mfilename, npairs, nunique_pairs);
end

%-----------------------------------------------------------------------------
% quick test
%-----------------------------------------------------------------------------
function test_aggregate_pairs
p = [[1, 2]; [1, 2]; [1, 5]; [2, 1]; [2, 6]; [1, 5]; [3, 4]; [1, 2]; [2, 6]];
pw = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9]; 

[op, osp, oi] = aggregate_pairs(p, pw);

p_check = op(oi, :);
pw_check = osp(oi);

assert(all(vec(p_check == p)), 'constraints are not the same');
assert(all(vec(pw_check >= pw)), 'problem with weights');

fprintf('%s: tests passed\n', mfilename);
end
