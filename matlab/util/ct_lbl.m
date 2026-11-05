function s = ct_lbl(s)
%CT_LBL  Axis/legend text: the degree sign is replaced in Octave (gnuplot cannot always render it).
if ct_is_octave()
    s = strrep(s, char(176), ' deg');
    s = strrep(s, char([194 176]), ' deg');
end
end
