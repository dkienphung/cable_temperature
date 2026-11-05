%RUN_PAPER_STUDY  Reproduce the study of Riba & Llaurado, Materials 2022, 15, 6814, in MATLAB/Octave.
%
%   >> run_paper_study                       % figures and summary in results/paper/
%
%   1. Test sequences of the paper (Tables 2 and 4) with Eq. 25 as printed and corrected,
%      compared with the digitized Figs. 9 and 11 (incl. relative difference panels)
%   2. Radial temperature profiles (Fig. 10a) and heat-balance terms (Fig. 10b)
%   3. Steady-state current-temperature relationship and ampacity vs IEC 60287
%   4. 1 s short circuit vs IEC 60949 (adiabatic and non-adiabatic)
%   5. Laboratory environment: elevation and indoor air movement
%   6. Diagnostic inverse fit (heat-capacity and h multipliers) with fminsearch
%
%   Same model and numbers as the Python package (cable_thermal); run_tests checks the agreement.

clear; clc;
root = setup_paths();
out = fullfile(root, 'results', 'paper');
if ~exist(out, 'dir'), mkdir(out); end
P = ct_palette();
D = ct_json('paper_figures_digitized.json');
deg = char(176);
log_lines = {};
say = @(varargin) fprintf(varargin{:});

%% 1. Test sequences (Figs. 9 and 11)
cables = {ct_cable_preset('cable70'), ct_cable_preset('cable150')};
figkeys = {'fig9_cable70', 'fig11_cable150'};
figname = {'fig09_cable70_temperature_rise.png', 'fig11_cable150_temperature_rise.png'};
titles = {'Conductor of 70 mm2: model (Table 1, corrected Eq. 25) vs paper Fig. 9', ...
          'Conductor of 150 mm2: model (Table 3, corrected Eq. 25) vs paper Fig. 11'};
probes = {'conductor_insulation', 'insulation_jacket', 'jacket_air'};
plabel = {'Conductor-insulation', 'Insulation-jacket', 'Jacket-air'};
pcol = {P.orange, P.blue, P.grey};
R = struct();
for ic = 1:2
    c = cables{ic};
    for conv = {'paper', 'ieee738'}
        m = ct_model(c, [], 'convection', conv{1});
        r = ct_simulate(m, c.test_steps);
        R.(c.name).(conv{1}) = r;
        s = ct_paper_compare(r, figkeys{ic});
        say('%s, Eq. 25 %-8s: end %.2f / %.2f / %.2f degC (paper %.2f / %.2f / %.2f), RMS %.2f K, energy error %.1e W/m\n', ...
            c.name, conv{1}, s.conductor_insulation.end_model, s.insulation_jacket.end_model, s.jacket_air.end_model, ...
            s.conductor_insulation.end_paper, s.insulation_jacket.end_paper, s.jacket_air.end_paper, s.all_rms, ct_energy_error(r));
        log_lines{end + 1} = sprintf('%s Eq25=%s end_model=%.3f/%.3f/%.3f end_paper=%.3f/%.3f/%.3f rms=%.3f mean=%+.3f', ...
            c.name, conv{1}, s.conductor_insulation.end_model, s.insulation_jacket.end_model, s.jacket_air.end_model, ...
            s.conductor_insulation.end_paper, s.insulation_jacket.end_paper, s.jacket_air.end_paper, s.all_rms, s.all_mean); %#ok<SAGROW>
    end
    r = R.(c.name).ieee738;
    s = ct_paper_compare(r, figkeys{ic});
    f = ct_figure(7.2, 7.6);
    subplot(2, 1, 1); hold on; box on; grid on;
    for k = 1:3
        plot(r.t, ct_probe(r, probes{k}), '-', 'color', pcol{k}, 'linewidth', 1.8, 'displayname', [plabel{k} ' (model)']);
        Pk = s.(probes{k}).paper;
        plot(Pk(:, 1), Pk(:, 2), '--', 'color', pcol{k}, 'linewidth', 1.1, 'displayname', [plabel{k} ' (paper, digitized)']);
    end
    edges = [0; cumsum(c.test_steps(:, 2))];
    for k = 1:size(c.test_steps, 1)
        tm = 0.5 * (edges(k) + edges(k + 1));
        text(tm, interp1(r.t, ct_probe(r, 'conductor_insulation'), tm) + 0.6, sprintf('%g A', c.test_steps(k, 1)), ...
             'horizontalalignment', 'center', 'fontweight', 'bold');
    end
    xlim([0 r.t(end)]); ylabel(ct_lbl(['Temperature [' deg 'C]'])); title(titles{ic});
    ct_legend('southeast'); set(gca, 'fontsize', 9);
    subplot(2, 1, 2); hold on; box on; grid on;
    for k = 1:3
        Pk = s.(probes{k}).paper;
        Tm = interp1(r.t, ct_probe(r, probes{k}), Pk(:, 1));
        plot(Pk(:, 1), 100 * (Pk(:, 2) - Tm) ./ Pk(:, 2), '-', 'color', pcol{k}, 'linewidth', 1.3, 'displayname', plabel{k});
    end
    plot([0 r.t(end)], [0 0], 'k-', 'handlevisibility', 'off'); xlim([0 r.t(end)]);
    xlabel('Time [s]'); ylabel('(paper - model)/paper [%]'); ct_legend('southwest'); set(gca, 'fontsize', 9);
    ct_save_fig(f, fullfile(out, figname{ic}));
end

%% 2. Fig. 10a radial profiles and Fig. 10b heat balance (70 mm2)
c = cables{1};
m = ct_model(c);
r = R.cable70.ieee738;
f = ct_figure(7, 5); hold on; box on; grid on;
NA = m.opt.nodes(1); NB = m.opt.nodes(2);
nodes = 1:m.n;
for t = 0:1000:8000
    [~, j] = min(abs(r.t - t));
    plot(nodes, r.T(j, :), 'k-', 'linewidth', 1.1);
    jl = NA + round(NB / 2);
    text(jl, r.T(j, jl) + 0.3, sprintf('t = %d s', t), 'fontsize', 8);
end
yl = ylim;
plot([NA + 1 NA + 1], yl, ':', 'color', P.ink2, 'handlevisibility', 'off'); plot([NA + NB + 1 NA + NB + 1], yl, ':', 'color', P.ink2, 'handlevisibility', 'off');
text(NA / 2, yl(2), 'CONDUCTOR', 'verticalalignment', 'top', 'horizontalalignment', 'center');
text(NA + 1 + NB / 2, yl(2), 'INSULATION', 'verticalalignment', 'top', 'horizontalalignment', 'center');
text(NA + NB + 1 + (m.n - NA - NB - 1) / 2, yl(2), 'JACKET', 'verticalalignment', 'top', 'horizontalalignment', 'center');
xlabel('Node position'); ylabel(ct_lbl(['Temperature [' deg 'C]']));
title('Temporal evolution of the radial temperature (cf. paper Fig. 10a)');
ct_save_fig(f, fullfile(out, 'fig10a_cable70_radial_profiles.png'));

f = ct_figure(7, 5); hold on; box on; grid on;
t = r.t(2:end);
plot(t, r.joule(2:end), 'k-', 'linewidth', 2.2, 'displayname', 'Joule heating I^2r (model)');
plot(t, r.p_conv(2:end) + r.p_rad(2:end), 'k--', 'linewidth', 2, 'displayname', 'P_c + P_r (model)');
plot(t, r.p_conv(2:end), '-', 'color', P.orange, 'linewidth', 1.8, 'displayname', 'Convection P_c (model)');
plot(t, r.p_rad(2:end), '-', 'color', P.blue, 'linewidth', 1.8, 'displayname', 'Radiation P_r (model)');
hb = D.fig10b_cable70_heat_terms;
plot(hb.Pc(:, 1), hb.Pc(:, 2), ':', 'color', P.orange, 'linewidth', 1.8, 'displayname', 'P_c (paper Fig. 10b)');
plot(hb.Pr(:, 1), hb.Pr(:, 2), ':', 'color', P.blue, 'linewidth', 1.8, 'displayname', 'P_r (paper Fig. 10b)');
xlim([0 r.t(end)]); xlabel('Time [s]'); ylabel('Per unit length heat terms [W/m]'); ct_legend('northwest');
title('Components of the heat balance, 70 mm2 (cf. paper Fig. 10b)');
ct_save_fig(f, fullfile(out, 'fig10b_cable70_heat_balance.png'));
log_lines{end + 1} = sprintf('cable70 heat balance at 8600 s: joule=%.3f pc=%.3f pr=%.3f W/m', r.joule(end), r.p_conv(end), r.p_rad(end));

%% 2b. Eq. 25 printed vs corrected
f = ct_figure(7.2, 4.6); hold on; box on; grid on;
Pk = D.fig9_cable70.conductor_insulation;
plot(Pk(:, 1), Pk(:, 2), 'k--', 'linewidth', 1.4, 'displayname', 'Paper Fig. 9 (digitized)');
plot(R.cable70.paper.t, ct_probe(R.cable70.paper, 'conductor_insulation'), '-', 'color', P.red, 'linewidth', 1.8, ...
     'displayname', 'Model, Eq. (25) as printed (\rho^{0.25})');
plot(r.t, ct_probe(r, 'conductor_insulation'), '-', 'color', P.orange, 'linewidth', 1.8, ...
     'displayname', 'Model, Eq. (25) corrected (\rho^{0.5}, IEEE 738)');
xlabel('Time [s]'); ylabel(ct_lbl(['Conductor-insulation temperature [' deg 'C]'])); ct_legend('northwest');
title('70 mm2: effect of the Eq. (25) correction');
ct_save_fig(f, fullfile(out, 'fig09b_cable70_formula_variants.png'));

%% 3. Steady-state current-temperature relationship and ampacity (IEC 60287 check)
f = ct_figure(11, 4.3);
for ic = 1:2
    c = cables{ic};
    m = ct_model(c);
    Ta = c.ambient_temperature;
    Imax = ct_ampacity(m, 90, Ta);
    I = linspace(0, 1.15 * Imax, 30);
    Tc = zeros(size(I)); Ts = zeros(size(I));
    for k = 1:numel(I)
        T = ct_steady_state(m, I(k), Ta); Tc(k) = T(1); Ts(k) = T(end);
    end
    kB = ct_mat_k(c.insulation_material, 75); kC = ct_mat_k(c.jacket_material, 75);
    Ti = linspace(Ta + 1, 100, 30); Iiec = zeros(size(Ti)); Iiec2 = Iiec;
    for k = 1:numel(Ti)
        q = ct_iec60287_free_air(c, Ti(k), Ta, kB, kC); Iiec(k) = q.I;
        q = ct_iec60287_free_air(c, Ti(k), Ta, 1 / 3.5, 1 / 5.0); Iiec2(k) = q.I;
    end
    q90 = ct_iec60287_free_air(c, 90, Ta, kB, kC); q90s = ct_iec60287_free_air(c, 90, Ta, 1 / 3.5, 1 / 5.0);
    say('%s: ampacity at 90 degC = %.1f A (IEC 60287 with paper k: %.1f A, with IEC rho_T: %.1f A)\n', c.name, Imax, q90.I, q90s.I);
    log_lines{end + 1} = sprintf('%s ampacity90=%.2f iec_paper_k=%.2f iec_rhoT=%.2f', c.name, Imax, q90.I, q90s.I); %#ok<SAGROW>
    subplot(1, 2, ic); hold on; box on; grid on;
    plot(I, Tc, '-', 'color', P.orange, 'linewidth', 1.8, 'displayname', 'Conductor (model)');
    plot(I, Ts, '-', 'color', P.grey, 'linewidth', 1.8, 'displayname', 'Jacket surface (model)');
    plot(Iiec, Ti, '--', 'color', P.violet, 'linewidth', 1.5, 'displayname', 'IEC 60287 (k from paper curves)');
    plot(Iiec2, Ti, ':', 'color', P.aqua, 'linewidth', 1.8, 'displayname', 'IEC 60287 (IEC \rho_T)');
    plot([0 I(end)], [90 90], '-', 'color', P.red, 'handlevisibility', 'off');
    text(Imax, 93, ct_lbl(sprintf(['%.0f A @ 90 ' char(176) 'C'], Imax)), 'horizontalalignment', 'right', 'fontsize', 8);
    xlabel('Current [A rms]'); ylabel(ct_lbl(['Temperature [' deg 'C]']));
    title(ct_lbl(sprintf(['%s: steady current-temperature, T_{air} = %g ' char(176) 'C'], c.name, Ta))); ct_legend('northwest');
end
ct_save_fig(f, fullfile(out, 'fig13_steady_state_current_temperature.png'));

%% 4. Short circuit (70 mm2) vs IEC 60949
c = cables{1}; m = ct_model(c);
Iload = ct_ampacity(m, 90, c.ambient_temperature);
T0 = ct_steady_state(m, Iload, c.ambient_temperature);
beta_p = 1 / c.alpha - c.T0;
cp_ref = ct_mat_cp(c.conductor_material, 75);
K_p = sqrt(c.conductor_mass_per_length * cp_ref / (ct_cable_resistance(c, c.T0) * c.alpha)) / (c.effective_cross_section * 1e6);
tsc = 1;
Isc = ct_iec60949('current', 250, tsc, 70, T0(1), K_p, beta_p);
prof = @(t) Isc * (t < tsc) + Iload * (t >= tsc);
rsc = ct_simulate(m, prof, [], 60, T0, 0.005, 4);
r10 = ct_simulate(m, prof, [], 60, T0, 10);
eps_na = 1 + 0.41 * sqrt(tsc / 70) + 0.12 * tsc / 70;
T_na = ct_iec60949('final', Isc / eps_na, tsc, 70, T0(1), K_p, beta_p);
tt = linspace(0, tsc, 101);
T_ad = ct_iec60949('final', Isc, tt, 70, T0(1), K_p, beta_p);
say('short circuit %.2f kA, %g s from %.1f degC: model peak %.1f degC, IEC 60949 adiabatic %.1f, non-adiabatic %.1f, dt=10 s peak %.1f\n', ...
    Isc / 1e3, tsc, T0(1), max(ct_probe(rsc, 'center')), T_ad(end), T_na, max(ct_probe(r10, 'conductor_insulation')));
log_lines{end + 1} = sprintf('short_circuit Isc=%.1f K=%.2f beta=%.2f model_peak=%.2f iec_adiabatic=%.2f iec_nonadiabatic=%.2f', ...
    Isc, K_p, beta_p, max(ct_probe(rsc, 'center')), T_ad(end), T_na);
f = ct_figure(7.5, 4.6); hold on; box on; grid on;
for k = 1:3
    plot(rsc.t, ct_probe(rsc, probes{k}), '-', 'color', pcol{k}, 'linewidth', 1.8, 'displayname', [plabel{k} ' (dt = 5 ms)']);
end
plot(tt, T_ad, '--', 'color', P.violet, 'linewidth', 1.6, 'displayname', 'IEC 60949 adiabatic');
plot(tsc, T_na, 'd', 'color', P.aqua, 'markerfacecolor', P.aqua, 'markersize', 8, 'displayname', 'IEC 60949 non-adiabatic');
plot(r10.t, ct_probe(r10, 'conductor_insulation'), 's', 'color', P.red, 'markerfacecolor', P.red, 'displayname', 'paper time step dt = 10 s');
xlim([0 20]); xlabel('Time [s]'); ylabel(ct_lbl(['Temperature [' deg 'C]'])); ct_legend('east');
title(ct_lbl(sprintf(['70 mm2: %.2f kA, %g s short circuit from %.0f ' char(176) 'C'], Isc / 1e3, tsc, T0(1))));
ct_save_fig(f, fullfile(out, 'fig14_short_circuit.png'));

%% 5. Laboratory environment: indoor air movement (IEEE 738 forced convection), H = 280 m
f = ct_figure(11, 4.4);
speeds = [0 0.1 0.2 0.3];
scol = {P.orange, P.violet, P.aqua, P.blue};
for ic = 1:2
    c = cables{ic};
    subplot(1, 2, ic); hold on; box on; grid on;
    Pk = D.(figkeys{ic}).conductor_insulation;
    plot(Pk(:, 1), Pk(:, 2), 'k--', 'linewidth', 1.4, 'displayname', 'Paper (digitized)');
    for k = 1:numel(speeds)
        r = ct_simulate(ct_model(c, [], 'air_speed', speeds(k)), c.test_steps);
        s = ct_paper_compare(r, figkeys{ic});
        plot(r.t, ct_probe(r, 'conductor_insulation'), '-', 'color', scol{k}, 'linewidth', 1.6, ...
             'displayname', sprintf('V = %.1f m/s (RMS %.2f K)', speeds(k), s.all_rms));
        log_lines{end + 1} = sprintf('%s air_speed=%.2f end=%.3f rms=%.3f', c.name, speeds(k), s.conductor_insulation.end_model, s.all_rms); %#ok<SAGROW>
    end
    xlabel('Time [s]'); ylabel(ct_lbl(['Conductor-insulation [' deg 'C]'])); ct_legend('northwest');
    title(sprintf('%s, H = %g m, indoor air movement', c.name, c.elevation));
end
ct_save_fig(f, fullfile(out, 'fig19_lab_environment_air_speed.png'));

%% 6. Diagnostic inverse fit (70 mm2): multipliers on heat capacity and h that reproduce Fig. 9
c = cables{1};
obj = @(p) ct_fit_rms(c, p, 'fig9_cable70');
[pfit, rmsfit] = fminsearch(obj, [1 1], optimset('TolX', 2e-3, 'TolFun', 1e-3, 'MaxFunEvals', 120));
say('diagnostic fit 70 mm2: heat capacity x%.2f, h x%.2f -> RMS %.2f K\n', abs(pfit(1)), abs(pfit(2)), rmsfit);
log_lines{end + 1} = sprintf('diagnostic_fit capacity_scale=%.3f h_scale=%.3f rms=%.3f', abs(pfit(1)), abs(pfit(2)), rmsfit);

fid = fopen(fullfile(out, 'summary.txt'), 'w');
fprintf(fid, '%s\n', log_lines{:});
fclose(fid);
say('figures and summary.txt written to %s\n', out);
