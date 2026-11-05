function tf = ct_gnuplot()
%CT_GNUPLOT  True when Octave plots with the gnuplot toolkit (limited font / symbol support).
tf = ct_is_octave() && strcmp(graphics_toolkit(), 'gnuplot');
end
