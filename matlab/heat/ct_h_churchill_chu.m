function [h, Ra, Pr] = ct_h_churchill_chu(Ts, Ta, D, H)
%CT_H_CHURCHILL_CHU  h = Nu k/D with Nu from Churchill & Chu and IEEE 738 air properties.
if nargin < 4 || isempty(H), H = 0; end
c = ct_const();
Tf = 0.5 * (Ts + Ta);
rho = ct_air_density(Ts, Ta, H);
[mu, k, cp] = ct_air_props(Tf);
nu = mu ./ rho;
alpha = k ./ (rho .* cp);
beta = 1 ./ (Tf + c.KELVIN);
Ra = c.G * beta .* abs(Ts - Ta) .* D.^3 ./ (nu .* alpha);
Pr = mu .* cp ./ k;
h = ct_nusselt_cc(Ra, Pr) .* k ./ D;
end
