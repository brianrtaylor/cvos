%-----------------------------------------------------------------------------
%
% @param: constraints : N x 2 matrix with node 1 and node 2
%
% @return: ind_pairs : the indices in the output of the corresponding i'th 
%   term of the input pairs / pair_weights set
%-----------------------------------------------------------------------------
function [unique_pairs, score_pairs, ind_pairs] = aggregate_pairs_fast( pairs, pair_weights)

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

[unique_pairs, score_pairs, ind_pairs] = ...
            aggregate_pairs_loop_mex(pairs, unique_pairs, pair_weights, sorti );
fprintf('%s: number pruned edges: %d --> %d\n', mfilename, npairs, nunique_pairs);
end

% testing function / script
function test_aggregate_pairs
p = [[1,2];[1,2];[1,5];[2,1];[2,6];[1,5];[3,4];[1,2];[2,6]];
pw = [0.1 ;  0.2;  0.3;  0.4;  0.5;  0.6;  0.7;  0.8;  0.9]; 

[op, osp, oi] = aggregate_pairs(p, pw);

p_check = op(oi, :);
pw_check = osp(oi);

assert(all(vec(p_check == p)), 'constraints are not the same');
assert(all(vec(pw_check >= pw)), 'problem with weights');

fprintf('%s: tests passed\n', mfilename);
end
