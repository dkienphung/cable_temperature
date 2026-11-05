function r = ct_cable_resistance(c, T)
%CT_CABLE_RESISTANCE  Paper Eq. (2): r(T) = r_T0 [1 + alpha (T - T0)], r_T0 = rho_e(T0)/S_eff [ohm/m].
r = c.resistivity_T0 / c.effective_cross_section * (1 + c.alpha * (T - c.T0));
end
