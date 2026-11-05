function s = ct_lbl(s)
%CT_LBL  Axis/legend text: the degree sign is replaced when Octave plots with gnuplot.
if ct_gnuplot()
    s = strrep(s, char(176), ' deg');
    s = strrep(s, char([194 176]), ' deg');
end
end
