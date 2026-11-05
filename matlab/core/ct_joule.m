function q = ct_joule(m, T, I)
%CT_JOULE  I_i^2 r_i(T_i) of every node [W/m]; sums to I^2 r(T) for uniform T (Eqs. 5, 8, 11).
q = m.opt.joule_scale * I^2 * ct_cable_resistance(m.cable, T) * m.r_scale .* m.joule_share;
end
