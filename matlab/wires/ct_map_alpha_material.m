function key = ct_map_alpha_material(name)
%CT_MAP_ALPHA_MATERIAL  Alpha Wire insulation name -> material key ('' = no thermal data).
%   Same mapping as data/insulation_materials.json ("alpha_wire_mapping").
names = {'PVC', 'PVC, Irradiated', 'PVC, Semi Rigid', 'XLPE', 'PTFE', 'PTFE TAPE', 'ETFE', 'FLUOROPLASTIC', ...
         'PVDF', 'Silicone Rubber', 'MPPE', 'MPPE-TPE', 'EPR Rubber', 'NYLON', 'Mica Glass Tape', 'PTFE Impregnated Fiberglass'};
keys = {'PVC', 'PVC', 'PVC', 'XLPE', 'PTFE', 'PTFE', 'ETFE', 'fluoropolymer_unspecified', ...
        'PVDF', 'silicone_rubber', 'mPPE', 'mPPE', 'EPR', 'PVC', '', ''};
i = find(strcmp(names, name), 1);
if isempty(i), error('ct_map_alpha_material:name', 'no mapping for %s', name); end
key = keys{i};
end
