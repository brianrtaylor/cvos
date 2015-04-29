%----------------------------------------------------------------------------%
% s_out = structsubset(s, fields)
%
% s_out is a struct containing the values listed in fields from sin
% 
% @param: s : input struct
% @param: fields : the fields we want, can be any number of fields
% @return: s_out : output struct with merged fields
%----------------------------------------------------------------------------%
function s_out = structsubset(s, varargin)
s_out = struct;
% for f = fieldnames(s2)';
for f = varargin;
  s_out = setfield(s_out, f{1}, getfield(s, f{1}));
end
end
