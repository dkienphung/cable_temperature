function m = ct_material(name, density)
%CT_MATERIAL  Temperature-dependent material digitized from the paper (Figs. 2-3).
%   m = ct_material('copper'|'xlpe'|'pvc', density)
%   Fields: name, density [kg/m3], cpT/cpV (degC, J/(kg K)), kT/kV (degC, W/(m K)).
%   Values are interpolated linearly and held constant outside 15/20-75 degC.
defaults = struct('copper', 8960, 'xlpe', 920, 'pvc', 1400);
if nargin < 2 || isempty(density), density = defaults.(name); end
d = ct_json('material_properties_digitized.json');
cp = d.([name '_cp']).data;
k = d.([name '_k']).data;
m = struct('name', name, 'density', density, 'cpT', cp(:, 1)', 'cpV', cp(:, 2)', 'kT', k(:, 1)', 'kV', k(:, 2)');
end
