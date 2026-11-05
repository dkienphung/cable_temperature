function h = ct_surface_h(m, Ts, Ta)
%CT_SURFACE_H  Convective coefficient at the outer surface [W/(m2 K)] (Eq. 25 / options).
if m.opt.convection_off || m.opt.adiabatic, h = 0; return; end
h = m.opt.h_scale * ct_h_surface(Ts, Ta, m.geom.D_out, m.elevation, m.opt.convection, m.opt.air_speed, m.opt.air_angle);
end
