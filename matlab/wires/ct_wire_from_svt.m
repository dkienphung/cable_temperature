function w = ct_wire_from_svt(row, meta, packing)
%CT_WIRE_FROM_SVT  General Cable Super Vu-Tron welding cable row -> wire record.
%   Conductor diameter is not published: estimated from the metal area with a packing factor.
if nargin < 3, packing = 0.80; end
area = row.strands * pi / 4 * ct_awg_diameter(row.strand_awg)^2;
d = sqrt(4 * area / (pi * packing));
od = row.od_in * 25.4e-3;
w = ct_wire_record('General Cable', ['Super Vu-Tron ' row.catalog], row.awg, d, area, 1.02 * 1.7241e-8 / area, ...
                   'computed: annealed Cu x 1.02 stranding factor', {'EPR'}, (od - d) / 2, od, meta.temp_max_C, ...
                   'rubber (modelled as EPR)', {sprintf('conductor diameter estimated with packing factor %.2f', packing)}, ...
                   row.current_rating_A);
end
