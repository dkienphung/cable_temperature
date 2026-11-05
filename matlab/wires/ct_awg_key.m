function k = ct_awg_key(awg)
%CT_AWG_KEY  Numeric sort key of an AWG string (larger = thinner; '1/0' = 0, '2/0' = -1).
s = strtrim(awg);
if numel(s) > 2 && strcmp(s(end - 1:end), '/0'), k = 1 - str2double(s(1:end - 2)); else, k = str2double(s); end
end
