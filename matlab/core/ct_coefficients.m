function [a, b, cc, d, C, q, h, prad] = ct_coefficients(m, T, I, Ta, dt)
%CT_COEFFICIENTS  TDMA coefficients a, b, c, d of Eqs. (7), (10), (13), (16), (19), (22), (27).
%   dt = Inf gives the steady-state equations (capacity terms dropped, radiation linearized).
G = ct_face_conductance(m, T);
C = ct_heat_capacity(m, T);
if isfinite(dt), cap = C / dt; else, cap = zeros(1, m.n); end
n = m.n;
a = zeros(1, n); cc = zeros(1, n);
a(2:n) = -G; cc(1:n - 1) = -G;
b = cap;
b(2:n) = b(2:n) + G;
b(1:n - 1) = b(1:n - 1) + G;
q = ct_joule(m, T, I);
d = cap .* T + q;
Ts = T(n);
h = ct_surface_h(m, Ts, Ta);
b(n) = b(n) + h * m.p_out;                 % last node, Eq. (27)
d(n) = d(n) + h * m.p_out * Ta;
[prad, s] = ct_surface_rad(m, Ts, Ta);
if strcmp(m.opt.radiation, 'linearized') || ~isfinite(dt)
    b(n) = b(n) + s;
    d(n) = d(n) - prad + s * Ts;
else
    d(n) = d(n) - prad;                    % explicit (lagged) radiation, Eq. (23)
end
end
