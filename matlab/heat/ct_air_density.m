function rho = ct_air_density(Ts, Ta, H, film_in_kelvin)
%CT_AIR_DENSITY  Paper Eq. (26) = IEEE 738-2012 Eq. (14a), temperatures in degC.
%   rho = (1.293 - 1.525e-4 H + 6.379e-9 H^2) / (1 + 0.00367 T_film)
%   film_in_kelvin = true reproduces the literal (wrong) reading of the paper's text.
if nargin < 3 || isempty(H), H = 0; end
if nargin < 4, film_in_kelvin = false; end
Tf = 0.5 * (Ts + Ta);
if film_in_kelvin, Tf = Tf + 273.15; end
rho = (1.293 - 1.525e-4 * H + 6.379e-9 * H.^2) ./ (1 + 0.00367 * Tf);
end
