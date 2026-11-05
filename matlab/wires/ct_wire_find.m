function w = ct_wire_find(part)
%CT_WIRE_FIND  Look up an Alpha Wire part number (e.g. '3051', '2840/7', 'PV0430') or a Habia type ('E 2219').
P = ct_alpha_parts();
for i = 1:numel(P)
    if strcmpi(P{i}.part_number, strtrim(part)), w = ct_wire_from_alpha(P{i}); return; end
end
refs = ct_json('reference_ratings.json');
rows = ct_as_cell(refs.habia_E_ptfe.rows);
for i = 1:numel(rows)
    if strcmpi(strrep(rows{i}.type, ' ', ''), strrep(part, ' ', ''))
        w = ct_wire_from_habia(rows{i}, refs.habia_E_ptfe); return;
    end
end
error('ct_wire_find:part', 'part %s not found', part);
end
