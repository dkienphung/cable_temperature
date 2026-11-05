function [mu, k, cp] = ct_air_props(Tf)
%CT_AIR_PROPS  Air properties at film temperature Tf [degC].
%   mu: IEEE 738-2012 Eq. (13a) [Pa s];  k: Eq. (15a) [W/(m K)];
%   cp: Incropera & DeWitt Table A.4 [J/(kg K)].
mu = 1.458e-6 * (Tf + 273).^1.5 ./ (Tf + 383.4);
k = 2.424e-2 + 7.477e-5 * Tf - 4.407e-9 * Tf.^2;
cp = interp1([250 300 350 400 450 500], [1006 1007 1009 1014 1021 1030], min(max(Tf + 273.15, 250), 500));
end
