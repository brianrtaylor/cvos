%-----------------------------------------------------------------------------
%
% @param: constraints : N x 2 matrix with node 1 and node 2
%
% @return: ind_pairs : the indices in the output of the corresponding i'th 
%   term of the input pairs / pair_weights set
%-----------------------------------------------------------------------------
function [unique_out, score_out, ind_out] = aggregate_list( ...
  in, in_weights, POR)

if isempty(in) && isempty(in_weights)
  unique_out = in;
  score_out = in_weights;
  ind_out = [];
  return
end

SAFE = false;
if ~exist('POR', 'var'); POR = true; end; % probabilistic or

% test function 
if nargin == 0; help aggregate_list; return; end;

n = size(in, 1);

if ~exist('in_weights', 'var'); in_weights = ones(n, 1); end;

[in, sorti] = sort(in);
in_weights = in_weights(sorti);

unique_out = unique(in);
nunique_out = size(unique_out, 1);
score_out = zeros(nunique_out, 1);
ind_out= zeros(n, 1);

k = 1; sk = 1;
while k <= n;
  % if pairs and unique_pairs match, add pair score to compressed unique_pairs
  % and increment to see next pair. Else, increment to next unique_pair
  if (in(k) == unique_out(sk));
    if ~POR;
      score_out(sk) = score_out(sk) + in_weights(k);
    else
      if SAFE;
        score_out(sk) = safepor(score_out(sk), in_weights(k));
      else
        score_out(sk) = score_out(sk) + in_weights(k) ...
          - score_out(sk) * in_weights(k);
      end
    end

    ind_out(sorti(k)) = sk;
    k = k + 1;
  else
    sk = sk + 1;
  end
end

nunique_out = size(unique_out, 1);
fprintf('%s: number pruned edges: %d --> %d\n', mfilename, n, nunique_out);
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
