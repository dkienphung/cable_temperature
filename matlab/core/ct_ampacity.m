function I = ct_ampacity(m, Tmax, Ta, location)
%CT_AMPACITY  Steady-state current that gives Tmax at a node (default: conductor centre).
%   The upper bracket is grown gradually: far above the rating the linear r(T) law can lead to
%   thermal runaway (no steady state), which is treated as 'too hot'.
if nargin < 2 || isempty(Tmax), Tmax = m.cable.max_continuous_temperature; end
if nargin < 3 || isempty(Ta), Ta = m.cable.ambient_temperature; end
if nargin < 4 || isempty(location), location = 'center'; end
idx = m.probes.(location);
f = @(x) excess(m, x, Ta, Tmax, idx);
lo = 0; hi = 50; fh = f(hi);
while fh < 0
    lo = hi; hi = 1.5 * hi; fh = f(hi);
    if hi > 1e6, error('ct_ampacity:bracket', 'no bracket found'); end
end
while ~isfinite(fh)            % runaway at hi: move hi back towards lo
    mid = 0.5 * (lo + hi); fm = f(mid);
    if isfinite(fm) && fm < 0, lo = mid; else, hi = mid; fh = fm; end
end
I = fzero(f, [lo hi], optimset('TolX', 1e-6));
end

function e = excess(m, I, Ta, Tmax, idx)
try
    T = ct_steady_state(m, I, Ta);
    e = T(idx) - Tmax;
catch
    e = Inf;
end
end
