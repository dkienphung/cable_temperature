function lg = ct_legend(loc, fs)
%CT_LEGEND  Legend with a small font at the given location.
if nargin < 2, fs = 7; end
lg = legend('location', loc);
set(lg, 'fontsize', fs);
end
