function h = ct_h_forced(Ts, Ta, D, V, H, angle_deg)
%CT_H_FORCED  IEEE 738-2012 forced convection (Eqs. 3a, 3b, 4a) expressed as h [W/(m2 K)].
%   q_c1 = K [1.01 + 1.35 Re^0.52] k_f dT,  q_c2 = K 0.754 Re^0.6 k_f dT,  h = max(q_c1,q_c2)/(pi D dT)
if nargin < 5 || isempty(H), H = 0; end
if nargin < 6 || isempty(angle_deg), angle_deg = 90; end
if V <= 0, h = 0 * Ts; return; end
Tf = 0.5 * (Ts + Ta);
[mu, k] = ct_air_props(Tf);
Re = D .* ct_air_density(Ts, Ta, H) * V ./ mu;
phi = angle_deg * pi / 180;
Ka = 1.194 - cos(phi) + 0.194 * cos(2 * phi) + 0.368 * sin(2 * phi);
h = Ka * k ./ (pi * D) .* max(1.01 + 1.35 * Re.^0.52, 0.754 * Re.^0.6);
end
