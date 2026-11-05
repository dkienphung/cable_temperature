function [q, slope] = ct_surface_rad(m, Ts, Ta)
%CT_SURFACE_RAD  Radiation per unit length eps sigma p (Ts^4 - Ta^4) [W/m] and its derivative.
if m.opt.radiation_off || m.opt.adiabatic, q = 0; slope = 0; return; end
c = ct_const();
f = m.opt.rad_scale * m.cable.emissivity * c.SIGMA * m.p_out;
q = f * ((Ts + c.KELVIN)^4 - (Ta + c.KELVIN)^4);
slope = 4 * f * (Ts + c.KELVIN)^3;
end
