function e = ct_energy_error(r)
%CT_ENERGY_ERROR  max |stored - (Joule - Pc - Pr)| [W/m]; must be ~1e-12 (discrete conservation).
e = max(abs(r.stored(2:end) - (r.joule(2:end) - r.p_conv(2:end) - r.p_rad(2:end))));
end
