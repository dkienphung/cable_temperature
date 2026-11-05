function c = ct_wire_to_cable(w, ambient, emissivity, elevation)
%CT_WIRE_TO_CABLE  Wire record -> cable struct for ct_model.
%   One insulation layer is represented as two equal sub-layers of the same material (the paper's
%   three-material scheme with B = C); insulation + jacket map to B and C.
if nargin < 2 || isempty(ambient), ambient = 30; end
if nargin < 3 || isempty(emissivity), emissivity = 0.85; end
if nargin < 4 || isempty(elevation), elevation = 0; end
if ~w.simulatable, error('ct_wire_to_cable:data', '%s: no thermal data for %s', w.part_number, w.insulation_label); end
if numel(w.layer_keys) == 1
    mB = ct_insulation_material(w.layer_keys{1}); mC = mB; th = [w.layer_t / 2, w.layer_t / 2];
elseif numel(w.layer_keys) == 2
    mB = ct_insulation_material(w.layer_keys{1}); mC = ct_insulation_material(w.layer_keys{2}); th = w.layer_t;
else
    error('ct_wire_to_cable:layers', 'at most two non-conductor layers are supported');
end
tmax = w.temp_max_C; if isnan(tmax), tmax = 90; end
c = ct_cable('name', w.part_number, 'designation', sprintf('%s %s %s AWG %s', w.source, w.part_number, w.awg, w.insulation_label), ...
    'conductor_diameter', w.conductor_diameter, 'insulation_thickness', th(1), 'jacket_thickness', th(2), ...
    'effective_cross_section', w.metal_area, 'resistivity_T0', w.r20 * w.metal_area, 'alpha', 0.00393, 'T0', 20, ...
    'conductor_mass_per_length', 8960 * w.metal_area, 'conductor_material', ct_material('copper'), ...
    'insulation_material', mB, 'jacket_material', mC, 'emissivity', emissivity, ...
    'ambient_temperature', ambient, 'elevation', elevation, 'max_continuous_temperature', tmax);
end
