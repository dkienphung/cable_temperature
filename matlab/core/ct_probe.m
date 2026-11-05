function v = ct_probe(r, name)
%CT_PROBE  Temperature history at a measuring point: 'center', 'conductor_insulation',
%   'insulation_jacket', 'jacket_air'.
v = r.T(:, r.probes.(name));
end
