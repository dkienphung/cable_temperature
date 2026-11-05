function g = ct_cable_geometry(c)
%CT_CABLE_GEOMETRY  Radii, outer diameter, conductor region area, mass and homogenized density.
g.R_A = 0.5 * c.conductor_diameter;
g.R_B = g.R_A + c.insulation_thickness;
g.R_C = g.R_B + c.jacket_thickness;
g.D_out = 2 * g.R_C;
g.S_geo = pi * g.R_A^2;
if isnan(c.conductor_mass_per_length)
    g.mass = c.conductor_material.density * c.effective_cross_section;
else
    g.mass = c.conductor_mass_per_length;
end
g.rho_eff = g.mass / g.S_geo;   % so that sum(rho_A S_i,A) = m  (Eq. 1)
end
