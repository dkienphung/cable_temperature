function w = ct_wire_from_rs(row, meta)
%CT_WIRE_FROM_RS  RS PRO PTFE type A row -> wire record (OD = mid-point of the published min/max).
d = row.conductor_diameter_mm * 1e-3;
od = 0.5 * (row.od_mm_min + row.od_mm_max) * 1e-3;
p = strsplit(row.stranding, '/');
area = str2double(p{1}) * pi / 4 * (str2double(p{2}) * 1e-3)^2;
w = ct_wire_record('RS PRO', row.stranding, row.awg, d, area, row.dcr_ohm_per_km_20C_max * 1e-3, ...
                   'datasheet DCR (maximum)', {'PTFE'}, (od - d) / 2, od, meta.temp_max_C, 'PTFE', {}, row.current_rating_A);
end
