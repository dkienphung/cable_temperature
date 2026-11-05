function f = ct_figure(w, h)
%CT_FIGURE  Invisible figure of w x h inches (headless-safe).  In Octave the Qt toolkit is used
%   when available (e.g. under xvfb-run), otherwise gnuplot.
if ct_is_octave()
    tk = available_graphics_toolkits();
    if any(strcmp(tk, 'qt')), graphics_toolkit('qt'); else, graphics_toolkit('gnuplot'); end
end
f = figure('visible', 'off', 'units', 'inches', 'position', [0 0 w h], 'paperpositionmode', 'auto', 'color', 'w');
end
