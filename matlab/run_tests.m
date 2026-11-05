%RUN_TESTS  Verification of the MATLAB/Octave implementation.
%
%   >> run_tests
%
%   Formula and solver checks (TDMA, mesh areas of Eqs. 6-24, energy conservation, analytic
%   steady state, adiabatic short circuit vs IEC 60949, Eq. 25 vs IEEE 738, Churchill-Chu,
%   IEC 60287 iteration) and parity with the Python package (test sequences of Figs. 9 and 11,
%   wire ampacities of the database).  Prints PASS/FAIL per test; errors at the end if any fail.

clear; clc;
root = setup_paths();
check = @(name, ok, msg) fprintf('%-62s %s  %s\n', name, ct_passfail(ok), msg);
nfail = 0;

%% TDMA vs direct solver (Eqs. 3-4)
n = 28; a = -rand(1, n); c = -rand(1, n); a(1) = 0; c(end) = 0; b = 3 + rand(1, n); d = rand(1, n);
A = diag(b) + diag(a(2:end), -1) + diag(c(1:end - 1), 1);
e = max(abs(ct_tdma(a, b, c, d) - (A \ d')'));
ok = e < 1e-12; nfail = nfail + ~ok; check('TDMA = dense solve', ok, sprintf('max diff %.1e', e));
e = max(abs(ct_solve_tridiag(a, b, c, d) - (A \ d')'));
ok = e < 1e-12; nfail = nfail + ~ok; check('sparse tridiagonal solve = dense solve', ok, sprintf('max diff %.1e', e));

%% Mesh of the paper (Eqs. 6, 9, 12, 15, 18, 21, 24)
for name = {'cable70', 'cable150'}
    cab = ct_cable_preset(name{1}); m = ct_model(cab); g = m.geom;
    ok = m.n == 28 && abs(sum(m.S(:, 1)) - pi * g.R_A^2) < 1e-15 && ...
         abs(sum(m.S(:, 2)) - pi * (g.R_B^2 - g.R_A^2)) < 1e-15 && abs(sum(m.S(:, 3)) - pi * (g.R_C^2 - g.R_B^2)) < 1e-15 && ...
         abs(m.S(1, 1) - pi * (0.5 * m.dx(1))^2) < 1e-18 && abs(sum(m.joule_share) - 1) < 1e-12 && ...
         abs(sum(m.rho(1) * m.S(:, 1)) - cab.conductor_mass_per_length) < 1e-12;
    nfail = nfail + ~ok; check(['mesh areas, perimeters, mass: ' name{1}], ok, sprintf('%d nodes', m.n));
end

%% Discrete energy conservation
cab = ct_cable_preset('cable70');
for rad = {'explicit', 'linearized'}
    r = ct_simulate(ct_model(cab, [], 'radiation', rad{1}), cab.test_steps);
    e = ct_energy_error(r); ok = e < 1e-8; nfail = nfail + ~ok;
    check(['energy balance (' rad{1} ' radiation)'], ok, sprintf('%.1e W/m', e));
end

%% Steady state vs analytic cylinder solution (constant properties)
cu = ct_material_constant('cu', 8960, 390, 398); ins = ct_material_constant('ins', 920, 2000, 0.33); jac = ct_material_constant('jac', 1400, 1000, 0.14);
c0 = cab; c0.conductor_material = cu; c0.insulation_material = ins; c0.jacket_material = jac; c0.alpha = 0;
for N = [9 40]
    m = ct_model(c0, [], 'nodes', [N N N]);
    T = ct_steady_state(m, 200);
    q = 200^2 * ct_cable_resistance(c0, 20); g = m.geom;
    TBC = T(end) + q / (2 * pi * 0.14) * log(g.R_C / g.R_B);
    TAB = TBC + q / (2 * pi * 0.33) * log(g.R_B / g.R_A);
    T0 = TAB + q / (4 * pi * 398);
    tol = [0.02 0.002]; tol = tol((N == 40) + 1) * (T0 - T(end));
    ok = abs(T(m.probes.insulation_jacket) - TBC) < tol && abs(T(m.probes.conductor_insulation) - TAB) < tol && abs(T(1) - T0) < tol;
    nfail = nfail + ~ok; check(sprintf('steady state vs analytic (N = %d)', N), ok, sprintf('dT_AB err %.1e K', abs(T(m.probes.conductor_insulation) - TAB)));
end

%% Adiabatic short circuit vs IEC 60949
rho20 = 1.7241e-8; beta = 234.5; sig = 3.45e6; S = 70e-6;
cs = cab; cs.resistivity_T0 = rho20; cs.alpha = 1 / (beta + 20); cs.effective_cross_section = S;
cs.conductor_mass_per_length = 8960 * S; cs.conductor_material = ct_material_constant('cu', 8960, sig / 8960, 398);
cs.insulation_material = ct_material_constant('i', 1e-9, 1, 1e-12); cs.jacket_material = cs.insulation_material;
K = ct_iec60949('K', sig, rho20, [], [], [], beta);
m = ct_model(cs, [], 'adiabatic', true);
r = ct_simulate(m, 10000, 90, 1, 90, 5e-4, 100);
Tf = ct_iec60949('final', 10000, 1, 70, 90, K, beta);
ok = abs(K - 226) < 0.5 && abs(r.T(end, 1) - Tf) < 0.05; nfail = nfail + ~ok;
check('adiabatic short circuit = IEC 60949', ok, sprintf('K = %.1f, model %.2f vs IEC %.2f degC', K, r.T(end, 1), Tf));

%% Surface heat transfer formulas
D = cab.conductor_diameter + 2 * (cab.insulation_thickness + cab.jacket_thickness); Ts = 45; Ta = 19;
qcn = 3.645 * ct_air_density(Ts, Ta)^0.5 * D^0.75 * (Ts - Ta)^1.25;           % IEEE 738-2012 Eq. (5c)
h = ct_h_natural(Ts, Ta, D, 0, 'ieee738');
ok = abs(h * pi * D * (Ts - Ta) - qcn) < 1e-10 && abs(ct_h_natural(Ts, Ta, D, 0, 'paper') * ct_air_density(Ts, Ta)^0.25 - h) < 1e-10;
nfail = nfail + ~ok; check('Eq. 25 corrected = IEEE 738 q_cn; printed differs by rho^0.25', ok, '');
ok = ct_air_density(30, 20, 0, true) < 0.7 && abs(ct_air_density(20, 20) - 1.293 / (1 + 0.00367 * 20)) < 1e-12;
nfail = nfail + ~ok; check('Eq. 26 needs degC (kelvin reading gives ~0.6 kg/m3)', ok, '');
nu = (0.6 + 0.387 * 1e4^(1/6) / (1 + (0.559 / 0.7)^(9/16))^(8/27))^2;
ok = abs(ct_nusselt_cc(1e4, 0.7) - nu) < 1e-12 && abs(nu - 4.366) < 2e-3 && abs(ct_nusselt_cc(0, 0.7) - 0.36) < 1e-12;
nfail = nfail + ~ok; check('Churchill-Chu Nusselt (hand value, Ra->0 limit)', ok, sprintf('Nu(1e4) = %.3f', nu));
ok = ct_h_churchill_chu(60, 20, 0.5e-3) > 1.6 * ct_h_natural(60, 20, 0.5e-3, 0, 'ieee738');
nfail = nfail + ~ok; check('thin wire: Churchill-Chu h > 1.6 x Eq. 25', ok, '');

%% IEC 60287 iteration self-consistency
q = ct_iec60287_free_air(cab, 90, 30, 1 / 3.5, 1 / 5.0);
W = q.I^2 * q.R;
ok = abs(W * q.T4 - q.dtheta_s) / q.dtheta_s < 1e-8 && abs(W * (q.T1 + q.T3 + q.T4) - 60) < 1e-6;
nfail = nfail + ~ok; check('IEC 60287 free-air iteration', ok, sprintf('I = %.1f A', q.I));

%% Parity with the Python package: test sequences (Figs. 9 and 11)
for name = {'cable70', 'cable150'}
    cab = ct_cable_preset(name{1});
    r = ct_simulate(ct_model(cab), cab.test_steps);
    Py = dlmread(fullfile(root, 'reference', ['python_' name{1} '_test_sequence.csv']), ',', 1, 0);
    e = max(max(abs(r.T(:, [1 10 19 28]) - Py(:, 4:7))));
    ok = e < 2e-4 && size(Py, 1) == numel(r.t); nfail = nfail + ~ok;
    check(['Python parity, test sequence ' name{1}], ok, sprintf('max |dT| %.1e K (CSV has 6 digits)', e));
end

%% Parity with the Python package: wire ampacities
fid = fopen(fullfile(root, 'reference', 'python_wire_results.csv'));
hdr = strsplit(fgetl(fid), ','); py = textscan(fid, repmat('%q', 1, numel(hdr)), 'Delimiter', ','); fclose(fid);
pcol = @(nm) py{find(strcmp(hdr, nm))};
pparts = pcol('part'); pI = str2double(pcol('I_rated_A')); pI2 = str2double(pcol('I_rated_A_eq25'));
test_parts = {'3051', '2840/7', '391000', '780001', 'PV0430', '6714', 'AZ221934', '2636', '5951', '1636'};
emax = 0;
for k = 1:numel(test_parts)
    w = ct_wire_find(test_parts{k}); c = ct_wire_to_cable(w, 30); j = find(strcmp(pparts, w.part_number));
    I1 = ct_ampacity(ct_model(c, [], 'convection', 'churchill_chu'), min(w.temp_max_C, 300), 30);
    I2 = ct_ampacity(ct_model(c, [], 'convection', 'ieee738'), min(w.temp_max_C, 300), 30);
    emax = max([emax, abs(I1 / pI(j) - 1), abs(I2 / pI2(j) - 1)]);
end
ok = emax < 1e-5; nfail = nfail + ~ok;
check(sprintf('Python parity, ampacity of %d database wires', numel(test_parts)), ok, sprintf('max rel. diff %.1e', emax));

%% Simulink step function (no Simulink needed): emulate the discrete loop of build_cable_simulink_model
cab = ct_cable_preset('cable70'); m = ct_model(cab);
assignin('base', 'CT_MODEL', m); clear ct_simulink_step
edges = [0; cumsum(cab.test_steps(:, 2))]; T = cab.ambient_temperature * ones(m.n, 1);
nsteps = round(edges(end) / m.opt.dt); Tlog = zeros(nsteps + 1, m.n); Tlog(1, :) = T';
for k = 1:nsteps
    t = (k - 1) * m.opt.dt;
    I = cab.test_steps(find(edges(1:end - 1) <= t + m.opt.dt / 2, 1, 'last'), 1);   % 1-D Lookup Table (flat)
    T = ct_simulink_step(T, I, cab.ambient_temperature);                             % MATLAB Function block
    Tlog(k + 1, :) = T';
end
r = ct_simulate(m, cab.test_steps);
e = max(abs(Tlog(:) - r.T(:))); ok = e < 1e-12; nfail = nfail + ~ok;
check('Simulink step (Unit Delay + lookup + ct_simulink_step) = ct_simulate', ok, sprintf('max |dT| %.1e K', e));

%% Database contents
P = ct_alpha_parts(); awgs = cellfun(@(p) p.awg, P, 'uniformoutput', false);
ok = numel(P) == 337 && all(ismember({'36', '30', '22', '10', '4', '1', '1/0'}, awgs)) && ...
     abs(ct_awg_diameter('1/0') - 8.252e-3) / 8.252e-3 < 1e-3 && abs(ct_awg_diameter(36) - 0.127e-3) < 1e-12;
[nn, dd] = ct_parse_stranding('133(7x19)/29', '8'); ok = ok && nn == 133 && abs(dd - ct_awg_diameter(29)) < 1e-15;
nfail = nfail + ~ok; check('wire database: 337 parts, 36 ... 1/0 AWG, AWG/stranding parsing', ok, '');

fprintf('\n%d test(s) failed\n', nfail);
if nfail > 0, error('run_tests:failed', '%d test(s) failed', nfail); end
