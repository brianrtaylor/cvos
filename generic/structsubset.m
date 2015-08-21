%-----------------------------------------------------------------------------
% s_out = structsubset(s, fields)
%
% s_out is a struct containing the values listed in fields from sin
% 
% @return: s_out: output struct with merged fields
% @param: s: input struct
% @param: fields: the fields we want, can be any number of fields
%-----------------------------------------------------------------------------
function s_out = structsubset(s, varargin)
s_out = struct;
for f = varargin;
  s_out = setfield(s_out, f{1}, getfield(s, f{1}));
end
end
