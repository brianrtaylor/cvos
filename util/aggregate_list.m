%-----------------------------------------------------------------------------
% aggregate_list
%
% given a list of pairs (or edges), this function will combine repeats into 
% a single entry and aggregate their associated weights depending on POR
%
% @return: unique_out: the remaining unique pairs
% @return: w_out: the weights associated with the remaining unique pairs
% @return: ind_out: indices of in corresponding to the elements of unique_out
% @param: in (Nx2): edges between node 1 (1st column) and node 2 (2nd column)
% @param: w_in (N): weights associated with each edge
% @param: POR: flag to combine weights via "probabilistic or" or "addition"
%-----------------------------------------------------------------------------
function [unique_out, w_out, ind_out] = aggregate_list(in, w_in, POR)

% parameters and input validation
if nargin == 0; test_aggregate_pairs; return; end;

if isempty(in) && isempty(w_in);
  unique_out = in;
  w_out = w_in;
  ind_out = [];
  return
end

SAFE = false;
n = size(in, 1);
if ~exist('POR', 'var'); POR = true; end; % probabilistic or
if ~exist('w_in', 'var'); w_in = ones(n, 1); end;

% do work
[in, sorti] = sort(in);
w_in = w_in(sorti);

unique_out = unique(in);
nunique_out = size(unique_out, 1);
w_out = zeros(nunique_out, 1);
ind_out= zeros(n, 1);

k = 1; sk = 1;
while k <= n;
  % if pairs and unique_pairs match, add pair score to compressed unique_pairs
  % and increment to see next pair. Else, increment to next unique_pair
  if (in(k) == unique_out(sk));
    if ~POR;
      w_out(sk) = w_out(sk) + w_in(k);
    else
      if SAFE;
        w_out(sk) = safepor(w_out(sk), w_in(k));
      else
        w_out(sk) = w_out(sk) + w_in(k) - w_out(sk) * w_in(k);
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

%-----------------------------------------------------------------------------
% testing script
%-----------------------------------------------------------------------------
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
