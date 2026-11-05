function ct_save_fig(fig, path)
%CT_SAVE_FIG  Save a figure as PNG in MATLAB and in Octave (size taken from the figure, inches).
d = fileparts(path);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
if ct_is_octave()
    pos = get(fig, 'position');                     % inches (ct_figure)
    print(fig, '-dpngcairo', sprintf('-S%d,%d', round(130 * pos(3)), round(130 * pos(4))), path);
else
    print(fig, '-dpng', '-r200', path);
end
close(fig);
end
