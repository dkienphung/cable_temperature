function w = ct_wire_from_habia(row, meta)
%CT_WIRE_FROM_HABIA  Habia type E PTFE row -> wire record.
d = row.conductor_diameter_mm * 1e-3; od = row.od_mm * 1e-3;
area = row.strands * pi / 4 * (row.strand_diameter_mm * 1e-3)^2;
w = ct_wire_record('Habia', row.type, row.awg, d, area, row.dcr_ohm_per_km_20C * 1e-3, 'datasheet DCR', ...
                   {'PTFE'}, (od - d) / 2, od, meta.temp_max_C, 'PTFE', {}, row.current_rating_A);
end
