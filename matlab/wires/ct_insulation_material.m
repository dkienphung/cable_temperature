function mat = ct_insulation_material(key)
%CT_INSULATION_MATERIAL  Material struct for a key of data/insulation_materials.json
%   (PVC and XLPE use the paper's temperature-dependent curves; the others are constant).
db = ct_json('insulation_materials.json');
p = db.materials.(key);
if strcmp(p.model, 'paper_curves')
    mat = ct_material(lower(key), p.density);
else
    mat = ct_material_constant(key, p.density, p.cp, p.k);
end
end
