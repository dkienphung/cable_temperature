function tau = ct_time_constant(m, I, Ta)
%CT_TIME_CONSTANT  Time for the conductor to reach 63.2 % of its steady rise after a current step.
Tss = ct_steady_state(m, I, Ta);
rise = Tss(1) - Ta;
C = sum(ct_heat_capacity(m, Ta * ones(1, m.n)));
Ts = Tss(end); D = m.geom.D_out;
G = ct_surface_h(m, Ts, Ta) * pi * D + ct_q_radiation(Ts, Ta, D, m.cable.emissivity) / max(Ts - Ta, 1e-9);
tau_est = C / G;
r = ct_simulate(m, I, Ta, 6 * tau_est, [], tau_est / 300);
y = r.T(:, 1) - Ta;
j = find(y >= 0.632 * rise, 1);
if isempty(j) || j == 1, tau = NaN; return; end
tau = interp1(y(j - 1:j), r.t(j - 1:j), 0.632 * rise);
end
