%----------------------------------------------------------------------------%
% structmerge(s1, s2)
%
% merges structures s1 and s2, replacing all fields in s1 with values in s2
% when a collision occurs. If subfields are structs, it's done recursively
% 
% @param: s1 : input struct
% @param: s2 : input struct whose fields win in collision
% @return: sout : output struct with merged fields
%----------------------------------------------------------------------------%
function sout = structmerge(s1, s2)
sout = s1;
for f = fieldnames(s2)';
  id = f{1};
  o2 = s2.(id);
  if isfield(s1, id);
    o1 = s1.(id);
    if isstruct(o1) && isstruct(o2);
      o2 = structmerge(o1, o2);
    end
  end
  sout.(f{1}) = o2;
end
end
