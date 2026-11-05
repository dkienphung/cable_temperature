function f = ct_figure(w, h)
%CT_FIGURE  Invisible figure of w x h inches (headless-safe).
if ct_is_octave()
    try, graphics_toolkit('gnuplot'); catch, end %#ok<CTCH>
end
f = figure('visible', 'off', 'units', 'inches', 'position', [0 0 w h], 'paperpositionmode', 'auto', 'color', 'w');
end
