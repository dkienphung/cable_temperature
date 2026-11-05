function [n, d] = ct_parse_stranding(s, awg)
%CT_PARSE_STRANDING  Number of strands and strand diameter [m] from '7/30', 'SOLID',
%   '133(7x19)/29', '259(37x7)25', '19/.0591', '665 (7x95)/30'.
s = upper(strrep(strtrim(s), ' ', ''));
if strcmp(s, 'SOLID'), n = 1; d = ct_awg_diameter(awg); return; end
tok = regexp(s, '^(\d+)(?:\([^)]*\))?/?(\.?\d+(?:\.\d+)?)$', 'tokens', 'once');
if isempty(tok), error('ct_parse_stranding:format', 'cannot parse stranding %s', s); end
n = str2double(tok{1});
if tok{2}(1) == '.', d = str2double(tok{2}) * 25.4e-3; else, d = ct_awg_diameter(str2double(tok{2})); end
end
