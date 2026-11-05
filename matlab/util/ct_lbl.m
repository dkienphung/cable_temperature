function s = ct_lbl(s)
%CT_LBL  Text with a degree sign written as char(176): unchanged in MATLAB; in Octave converted to
%   UTF-8 (Qt toolkit) or to ' deg' (gnuplot toolkit, which cannot always render it).
if ~ct_is_octave(), return; end
s = strrep(s, char([194 176]), char(176));
if ct_gnuplot()
    s = strrep(s, char(176), ' deg');
else
    s = strrep(s, char(176), char([194 176]));
end
end
