function d = ct_awg_diameter(awg)
%CT_AWG_DIAMETER  ASTM B258: d = 0.127 mm x 92^((36 - n)/39); '1/0' -> n = 0, '2/0' -> -1, ...
if ischar(awg)
    s = upper(strtrim(awg));
    if numel(s) > 2 && strcmp(s(end - 1:end), '/0')
        n = 1 - str2double(s(1:end - 2));
    else
        n = str2double(s);
    end
else
    n = awg;
end
d = 0.127e-3 * 92.^((36 - n) / 39);
end
