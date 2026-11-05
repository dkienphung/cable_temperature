function res = ct_iec60287_free_air(c, theta_c, theta_amb, k_ins, k_jac, Z, E, gexp)
%CT_IEC60287_FREE_AIR  IEC 60287 rating of a single LV cable in free air (no sun, no sheath).
%   I = sqrt(dtheta / (R (T1 + T3 + T4))),  T1 = rho_ins/(2 pi) ln(1 + 2 t1/dc),
%   T3 = rho_jac/(2 pi) ln(1 + 2 t3/Da),    T4 = 1/(pi De h dtheta_s^(1/4)),  h = Z/De^g + E
%   (single cable: Z = 0.21, E = 3.94, g = 0.60); dtheta_s from the IEC iteration.
if nargin < 6, Z = 0.21; E = 3.94; gexp = 0.60; end
geo = ct_cable_geometry(c);
De = geo.D_out;
h = Z / De^gexp + E;
T1 = 1 / k_ins / (2 * pi) * log(1 + 2 * c.insulation_thickness / c.conductor_diameter);
T3 = 1 / k_jac / (2 * pi) * log(1 + 2 * c.jacket_thickness / (2 * geo.R_B));
dtheta = theta_c - theta_amb;
KA = pi * De * h * (T1 + T3);
s = 2;
for it = 1:200
    sn = (dtheta / (1 + KA * s))^0.25;
    if abs(sn - s) < 1e-12, s = sn; break; end
    s = sn;
end
T4 = 1 / (pi * De * h * s);
R = ct_cable_resistance(c, theta_c);
res = struct('I', sqrt(dtheta / (R * (T1 + T3 + T4))), 'T1', T1, 'T3', T3, 'T4', T4, 'h', h, ...
             'dtheta_s', s^4, 'R', R);
end
