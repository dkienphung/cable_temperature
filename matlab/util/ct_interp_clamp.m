function v = ct_interp_clamp(x, y, T)
%CT_INTERP_CLAMP  Piecewise-linear interpolation, value held constant outside the table.
if numel(x) == 1
    v = y(1) * ones(size(T));
    return
end
T = min(max(T, x(1)), x(end));
v = interp1(x, y, T, 'linear');
end
