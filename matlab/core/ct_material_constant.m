function m = ct_material_constant(name, density, cp, k)
%CT_MATERIAL_CONSTANT  Constant-property material.
m = struct('name', name, 'density', density, 'cpT', 0, 'cpV', cp, 'kT', 0, 'kV', k);
end
