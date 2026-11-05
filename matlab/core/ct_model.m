function m = ct_model(cable, opt, varargin)
%CT_MODEL  Build the radial finite-difference model of the paper (Figs. 4-7, Eqs. 5-24).
%   m = ct_model(cable)                 default options (ct_options)
%   m = ct_model(cable, opt)            options struct
%   m = ct_model(cable, [], 'convection', 'churchill_chu', ...)  name/value overrides
%
%   Node 1 = conductor centre, node NA+1 = conductor-insulation boundary, node NA+NB+1 =
%   insulation-jacket boundary, node NA+NB+NC+1 = jacket-air boundary.  Each node owns the half
%   elements on both sides (areas S_i,A, S_i,B, S_i,C) and the faces have perimeter
%   p = 2 pi (x_i +- 0.5 dx), exactly as in Eqs. (6), (9), (12), (15), (18), (21), (24).
if nargin < 2 || isempty(opt), opt = ct_options(); end
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i + 1}; end
g = ct_cable_geometry(cable);
NA = opt.nodes(1); NB = opt.nodes(2); NC = opt.nodes(3);
dxA = g.R_A / NA; dxB = cable.insulation_thickness / NB; dxC = cable.jacket_thickness / NC;
x = [(0:NA) * dxA, g.R_A + (1:NB) * dxB, g.R_B + (1:NC) * dxC];
x(NA + 1) = g.R_A; x(NA + NB + 1) = g.R_B; x(end) = g.R_C;
n = numel(x);
seg_mat = [ones(1, NA), 2 * ones(1, NB), 3 * ones(1, NC)];
seg_dx = [dxA * ones(1, NA), dxB * ones(1, NB), dxC * ones(1, NC)];
S = zeros(n, 3);
for s = 1:n - 1
    xm = x(s) + 0.5 * seg_dx(s);
    S(s, seg_mat(s)) = S(s, seg_mat(s)) + pi * (xm^2 - x(s)^2);              % right half of node s
    S(s + 1, seg_mat(s)) = S(s + 1, seg_mat(s)) + pi * (x(s + 1)^2 - xm^2);  % left half of node s+1
end
m.cable = cable; m.opt = opt; m.geom = g;
m.x = x; m.n = n; m.seg_mat = seg_mat; m.seg_dx = seg_dx;
m.face_p = 2 * pi * (x(1:end - 1) + 0.5 * seg_dx);
m.S = S; m.dx = [dxA dxB dxC];
m.probes = struct('center', 1, 'conductor_insulation', NA + 1, 'insulation_jacket', NA + NB + 1, 'jacket_air', n);
m.joule_share = S(:, 1)' / g.S_geo;
m.rho = [g.rho_eff, cable.insulation_material.density, cable.jacket_material.density];
m.mats = {cable.conductor_material, cable.insulation_material, cable.jacket_material};
switch opt.joule_area
    case 'effective', m.r_scale = 1;
    case 'geometric', m.r_scale = cable.effective_cross_section / g.S_geo;
    otherwise, error('ct_model:joule_area', 'joule_area must be effective or geometric');
end
if isempty(opt.elevation), m.elevation = cable.elevation; else, m.elevation = opt.elevation; end
m.p_out = 2 * pi * g.R_C;
end
