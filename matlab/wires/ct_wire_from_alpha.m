function w = ct_wire_from_alpha(rec)
%CT_WIRE_FROM_ALPHA  Alpha Wire datasheet record -> wire record (see docs/WIRE_DATABASE.md).
IN = 25.4e-3; KFT = 304.8;
[n, ds] = ct_parse_stranding(rec.stranding, rec.awg);
area = n * pi / 4 * ds^2;
d = rec.conductor_diameter_in * IN;
notes = {};
if isfield(rec, 'dcr_ohm_per_kft_20C') && ~isempty(rec.dcr_ohm_per_kft_20C)
    r20 = rec.dcr_ohm_per_kft_20C / KFT; rsrc = 'datasheet DCR';
else
    r20 = 1.7241e-8 / area; rsrc = 'computed (no DCR on datasheet)';
    notes{end + 1} = 'DC resistance computed from area with annealed-copper resistivity';
end
L = ct_as_cell(rec.insulation_layers);
keys = cell(1, numel(L)); t = zeros(1, numel(L)); names = cell(1, numel(L));
for i = 1:numel(L)
    keys{i} = ct_map_alpha_material(L{i}.material); t(i) = L{i}.wall_in * IN; names{i} = L{i}.material;
end
label = strjoin(names, ' / ');
od = NaN;
if isfield(rec, 'insulation_od_in') && ~isempty(rec.insulation_od_in), od = rec.insulation_od_in; end
if isfield(rec, 'jacket') && ~isempty(rec.jacket)
    keys{end + 1} = ct_map_alpha_material(rec.jacket.material); t(end + 1) = rec.jacket.wall_in * IN;
    label = [label ' + ' rec.jacket.material ' jacket'];
    od = rec.jacket.od_in;
end
if isnan(od), od = d + 2 * sum(t); else, od = od * IN; end
[t, dev] = ct_fit_layers(d, od, t);
if dev > 0.05, notes{end + 1} = sprintf('outer wall adjusted by %.0f%% to match the datasheet OD', 100 * dev); end
tmax = NaN; if isfield(rec, 'temp_max_C'), tmax = rec.temp_max_C; end
rating = NaN; if isfield(rec, 'current_rating_A'), rating = rec.current_rating_A; end
w = ct_wire_record('Alpha Wire', rec.part_number, rec.awg, d, area, r20, rsrc, keys, t, od, tmax, label, notes, rating);
w.family = ''; if isfield(rec, 'family') && ~isempty(rec.family), w.family = rec.family; end
end
