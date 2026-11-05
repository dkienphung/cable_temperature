function [t, dev] = ct_fit_layers(d, od, t)
%CT_FIT_LAYERS  Adjust the outermost wall so that d + 2 sum(t) equals the datasheet OD.
inner = sum(t(1:end - 1));
t_last = (od - d) / 2 - inner;
if t_last <= 0
    dev = abs(d + 2 * sum(t) - od) / od;
    return
end
dev = abs(t_last - t(end)) / max(t(end), 1e-12);
t(end) = t_last;
end
