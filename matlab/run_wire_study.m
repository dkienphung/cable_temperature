%RUN_WIRE_STUDY  Apply the cable model to the wire database (Alpha Wire, 36 AWG ... 1/0 AWG) and
%   compare with manufacturer current ratings (Habia, RS PRO, Alpha PV wire, Super Vu-Tron).
%
%   >> run_wire_study                 % results/wires/ : CSV, comparisons.txt, figures
%
%   Settings below: ambient, convection law, time constants on/off, part subset.
%   Model choices (docs/WIRE_DATABASE.md): single insulation = two equal sub-layers, outer wall
%   fitted to the datasheet OD, alpha = 0.00393 1/K, still air, emissivity 0.85.

clear; clc;
root = setup_paths();
out = fullfile(root, 'results', 'wires');
if ~exist(out, 'dir'), mkdir(out); end
P = ct_palette(); deg = char(176);

AMBIENT = 30;                    % degC
CONV = 'churchill_chu';          % natural convection for thin wires ('ieee738' = paper Eq. 25)
DO_TAU = ~ct_is_octave();        % thermal time constants (slow in Octave)
SUBSET = [];                     % [] = all parts, or e.g. {'3051','391000'}
REUSE = true;                    % reuse results/wires/campaign.mat if it exists (delete it to recompute)

%% 1. Campaign over the whole database
parts = ct_alpha_parts();
hdr = {'part', 'family', 'awg', 'insulation', 'group', 'T_max_C', 'D_cond_mm', 'OD_mm', 'wall_mm', ...
       'metal_area_mm2', 'R20_ohm_per_km', 'simulated', 'I_rated_A', 'I_rated_A_eq25', 'I_105C_A', ...
       'I_rise40K_A', 'tau_s', 'notes'};
rows = cell(numel(parts), numel(hdr));
W = cell(1, numel(parts)); res = nan(numel(parts), 5);
t0 = tic;
cache = fullfile(out, 'campaign.mat');
if REUSE && exist(cache, 'file') && isempty(SUBSET)
    S = load(cache); rows = S.rows; res = S.res; W = S.W;
    fprintf('campaign results loaded from %s\n', cache);
    parts = {};
end
for i = 1:numel(parts)
    w = ct_wire_from_alpha(parts{i}); W{i} = w;
    if ~isempty(SUBSET) && ~any(strcmpi(SUBSET, w.part_number)), continue; end
    g = ct_wire_group(w);
    v = nan(1, 5);
    if w.simulatable && ~isnan(w.temp_max_C)
        c = ct_wire_to_cable(w, AMBIENT);
        m = ct_model(c, [], 'convection', CONV);
        Tmax = min(w.temp_max_C, 300);
        v(1) = ct_ampacity(m, Tmax, AMBIENT);
        v(2) = ct_ampacity(ct_model(c, [], 'convection', 'ieee738'), Tmax, AMBIENT);
        v(3) = ct_ampacity(m, 105, AMBIENT);
        v(4) = ct_ampacity(m, AMBIENT + 40, AMBIENT);
        if DO_TAU, v(5) = ct_time_constant(m, v(4), AMBIENT); end
    end
    res(i, :) = v;
    rows(i, :) = {w.part_number, w.family, w.awg, w.insulation_label, g, w.temp_max_C, w.conductor_diameter * 1e3, ...
                  w.outer_diameter * 1e3, (w.outer_diameter - w.conductor_diameter) / 2 * 1e3, w.metal_area * 1e6, ...
                  w.r20 * 1e3, w.simulatable, v(1), v(2), v(3), v(4), v(5), strjoin(w.notes, '; ')};
    if mod(i, 25) == 0, fprintf('%d / %d parts (%.0f s)\n', i, numel(parts), toc(t0)); end
end
keep = ~cellfun(@isempty, rows(:, 1));
if ~isempty(parts) && isempty(SUBSET), save('-v7', cache, 'rows', 'res', 'W'); end
ct_write_csv(fullfile(out, 'wire_results_matlab.csv'), hdr, rows(keep, :));
fprintf('campaign done: %d parts in %.0f s\n', sum(keep), toc(t0));

%% 2. Ampacity vs AWG coloured by insulation group (rated temperature, 30 degC still air)
groups = {'PVC', 'XLPE (+PVC jacket)', 'Silicone rubber', 'PTFE', 'ETFE', 'Fluoroplastic (ThermoThin)', 'PVDF', 'mPPE (EcoWire)', 'EPR'};
awgs = unique(rows(keep, 3));
[~, o] = sort(-cellfun(@ct_awg_key, awgs)); awgs = awgs(o);
f = ct_figure(11, 5.5); hold on; box on; grid on;
for gi = 1:numel(groups)
    sel = find(keep & strcmp(rows(:, 5), groups{gi}) & ~isnan(res(:, 1)));
    if isempty(sel), continue; end
    [~, col] = ct_wire_group(W{sel(1)});
    xs = cellfun(@(a) find(strcmp(awgs, a)), rows(sel, 3)) + (gi - 5) * 0.06;
    semilogy(xs, res(sel, 1), 'o', 'color', col, 'markerfacecolor', col, 'markersize', 4, 'displayname', groups{gi});
end
set(gca, 'yscale', 'log', 'xtick', 1:numel(awgs), 'xticklabel', awgs, 'fontsize', 8);
xlabel('Conductor size [AWG]'); ylabel('Model ampacity [A]'); ct_legend('northwest');
title(ct_lbl(sprintf('Ampacity at the rated conductor temperature, still air %g %sC (%s)', AMBIENT, deg, CONV)));
ct_save_fig(f, fullfile(out, 'fig21_wire_ampacity_by_material.png'));

%% 3. Comparison with manufacturer current ratings
refs = ct_json('reference_ratings.json');
L = {};
% --- Habia PTFE type E, 'Amps @ 40 degC'
meta = refs.habia_E_ptfe; hrows = ct_as_cell(meta.rows);
H = nan(numel(hrows), 5);
for i = 1:numel(hrows)
    w = ct_wire_from_habia(hrows{i}, meta); c = ct_wire_to_cable(w, 40);
    mc = ct_model(c, [], 'convection', CONV); me = ct_model(c, [], 'convection', 'ieee738');
    Tc = ct_steady_state(mc, hrows{i}.current_rating_A, 40); Te = ct_steady_state(me, hrows{i}.current_rating_A, 40);
    H(i, :) = [hrows{i}.current_rating_A, Tc(1), Te(1), ct_ampacity(mc, 200, 40), ct_ampacity(me, 200, 40)];
    L{end + 1} = sprintf('Habia %-7s %2s AWG rating %2d A: T at rating %.0f (Eq.25 %.0f) degC, I(200 degC) %.1f (Eq.25 %.1f) A', ...
        hrows{i}.type, hrows{i}.awg, hrows{i}.current_rating_A, Tc(1), Te(1), H(i, 4), H(i, 5)); %#ok<SAGROW>
end
% ambient correction factors (200 degC conductor) for three sizes vs Habia's table
fn = fieldnames(meta.ambient_correction);
Tamb = cellfun(@(s) str2double(regexprep(s, '\D', '')), fn)'; [Tamb, o] = sort(Tamb); fn = fn(o);
pub = cellfun(@(s) meta.ambient_correction.(s), fn)';
der = nan(3, numel(Tamb)); types = {'E 3007', 'E 2219', 'E 1419'};
for k = 1:3
    w = ct_wire_find(types{k}); mc = ct_model(ct_wire_to_cable(w, 40), [], 'convection', CONV);
    I40 = ct_ampacity(mc, 200, 40);
    for j = 1:numel(Tamb), der(k, j) = ct_ampacity(mc, 200, Tamb(j)) / I40; end
end
L{end + 1} = ['Habia ambient factors (published): ' sprintf('%g:%.2f ', [Tamb; pub])];
for k = 1:3, L{end + 1} = ['model ' types{k} ': ' sprintf('%g:%.2f ', [Tamb; der(k, :)])]; end %#ok<SAGROW>
% --- RS PRO PTFE type A, rating = 40 K rise (20 degC air assumed)
meta = refs.rs_ptfe_type_a; rrows = ct_as_cell(meta.rows);
RS = nan(numel(rrows), 5);
for i = 1:numel(rrows)
    w = ct_wire_from_rs(rrows{i}, meta); c = ct_wire_to_cable(w, 20);
    mc = ct_model(c, [], 'convection', CONV); me = ct_model(c, [], 'convection', 'ieee738');
    Tc = ct_steady_state(mc, rrows{i}.current_rating_A, 20); Te = ct_steady_state(me, rrows{i}.current_rating_A, 20);
    RS(i, :) = [rrows{i}.current_rating_A, Tc(1) - 20, Te(1) - 20, ct_ampacity(mc, 60, 20), ct_ampacity(me, 60, 20)];
    L{end + 1} = sprintf('RS PRO %s AWG %-7s rating %g A: rise %.1f K (Eq.25 %.1f K), I(40 K) %.2f A', ...
        rrows{i}.awg, rrows{i}.stranding, RS(i, 1), RS(i, 2), RS(i, 3), RS(i, 4)); %#ok<SAGROW>
end
% --- Alpha Wire PV wire (XLPE + PVC jacket, 90 degC), datasheet current @ 30 degC
pv = {'PV0830', 'PV0630', 'PV0430'}; PV = nan(3, 4);
for i = 1:3
    w = ct_wire_find(pv{i}); c = ct_wire_to_cable(w, 30);
    mc = ct_model(c, [], 'convection', CONV); me = ct_model(c, [], 'convection', 'ieee738');
    Tc = ct_steady_state(mc, w.rating_A, 30);
    PV(i, :) = [w.rating_A, Tc(1), ct_ampacity(mc, 90, 30), ct_ampacity(me, 90, 30)];
    L{end + 1} = sprintf('Alpha %s %s AWG rating %d A: T at rating %.1f degC, I(90 degC) %.1f A (Eq.25 %.1f A)', ...
        pv{i}, w.awg, w.rating_A, PV(i, 2), PV(i, 3), PV(i, 4)); %#ok<SAGROW>
end
% --- General Cable Super Vu-Tron welding cable (90 degC, 30 degC): 6 AWG ... 4/0 AWG
meta = refs.general_cable_super_vutron; srows = ct_as_cell(meta.rows); SV = nan(numel(srows), 5);
for i = 1:numel(srows)
    I90 = zeros(1, 3); packs = [0.75 0.80 0.85];
    for k = 1:3
        mc = ct_model(ct_wire_to_cable(ct_wire_from_svt(srows{i}, meta, packs(k)), 30), [], 'convection', CONV);
        I90(k) = ct_ampacity(mc, 90, 30);
        if k == 2, Tc = ct_steady_state(mc, srows{i}.current_rating_A, 30); end
    end
    SV(i, :) = [srows{i}.current_rating_A, Tc(1), I90];
    L{end + 1} = sprintf('Super Vu-Tron %-3s AWG rating %d A: T at rating %.1f degC, I(90 degC) %.0f / %.0f / %.0f A (packing 0.75/0.80/0.85)', ...
        srows{i}.awg, SV(i, 1), SV(i, 2), I90); %#ok<SAGROW>
end
fid = fopen(fullfile(out, 'comparisons.txt'), 'w'); fprintf(fid, '%s\n', L{:}); fclose(fid);
fprintf('%s\n', L{:});

f = ct_figure(12, 5);
subplot(1, 2, 1); hold on; box on; grid on;
loglog(H(:, 1), H(:, 4), 'o', 'color', P.violet, 'markerfacecolor', P.violet, 'displayname', ct_lbl(['Habia PTFE (200 ' deg 'C @ 40 ' deg 'C)']));
loglog(RS(:, 1), RS(:, 4), 's', 'color', P.aqua, 'markerfacecolor', P.aqua, 'displayname', 'RS PRO PTFE (40 K rise)');
loglog(PV(:, 1), PV(:, 3), 'd', 'color', P.yellow, 'markerfacecolor', P.yellow, 'displayname', ct_lbl(['Alpha PV (90 ' deg 'C @ 30 ' deg 'C)']));
loglog(SV(:, 1), SV(:, 4), '^', 'color', P.grey, 'markerfacecolor', P.grey, 'displayname', ct_lbl(['Super Vu-Tron (90 ' deg 'C @ 30 ' deg 'C)']));
loglog([1.5 600], [1.5 600], 'k-', 'displayname', 'model = rating');
set(gca, 'xscale', 'log', 'yscale', 'log'); xlim([1.5 600]); ylim([1.5 600]);
xlabel('Manufacturer current rating [A]'); ylabel('Model current, same criterion [A]'); ct_legend('northwest');
title('Model (Churchill-Chu) vs manufacturer ratings');
subplot(1, 2, 2); hold on; box on; grid on;
plot(1:size(H, 1), H(:, 2), 'o-', 'color', P.violet, 'markerfacecolor', P.violet, 'displayname', 'Churchill-Chu');
plot(1:size(H, 1), H(:, 3), 'o--', 'color', P.violet, 'displayname', 'Eq. 25 (paper)');
plot([1 size(H, 1)], [200 200], '-', 'color', P.red, 'displayname', ct_lbl(['PTFE rating 200 ' deg 'C']));
set(gca, 'xtick', 1:size(H, 1), 'xticklabel', cellfun(@(r) strrep(r.type, 'E ', ''), hrows, 'uniformoutput', false), 'fontsize', 7);
xlabel('Habia type E'); ylabel(ct_lbl(['Model conductor temperature at rating [' deg 'C]'])); ct_legend('southeast');
title(ct_lbl(['Habia PTFE: temperature at the published current, 40 ' deg 'C air']));
ct_save_fig(f, fullfile(out, 'fig23_model_vs_manufacturer_ratings.png'));

f = ct_figure(7, 4.4); hold on; box on; grid on;
plot(Tamb, pub, 's-', 'color', P.ink, 'markerfacecolor', P.ink, 'displayname', 'Habia published');
dcol = {P.violet, P.aqua, P.orange};
for k = 1:3, plot(Tamb, der(k, :), 'o--', 'color', dcol{k}, 'displayname', ['model ' types{k}]); end
xlabel(ct_lbl(['Ambient air temperature [' deg 'C]'])); ylabel(ct_lbl(['Current correction factor (40 ' deg 'C = 1)'])); ct_legend('southwest');
title(ct_lbl(['Ambient correction: model (200 ' deg 'C conductor) vs Habia']));
ct_save_fig(f, fullfile(out, 'fig24_habia_ambient_correction.png'));

%% 4. Same 22 AWG conductor, different insulation: 8 A step in 30 degC still air
sel22 = {'3051', 'PVC 0.016" (UL1007)'; '3071', 'PVC 0.032" (UL1015)'; '5855/7', 'PTFE 0.010"'; ...
         '39X2205', 'Silicone 0.020"'; '6713', 'mPPE 0.0095"'; '7024', 'XLPE 0.016"'; ...
         'AZ221934', 'ETFE 0.011"'; '2622', 'Fluoroplastic 0.007"'};
cols = {P.orange, P.red, P.violet, P.blue, P.green, P.yellow, P.aqua, P.magenta};
f = ct_figure(8, 4.8); hold on; box on; grid on;
L22 = {};
for k = 1:size(sel22, 1)
    w = ct_wire_find(sel22{k, 1}); m = ct_model(ct_wire_to_cable(w, AMBIENT), [], 'convection', CONV);
    r = ct_simulate(m, 8, AMBIENT, 600, [], 0.5);
    plot(r.t, r.T(:, 1), '-', 'color', cols{k}, 'linewidth', 1.6, 'displayname', sprintf('%s %s, OD %.2f mm', sel22{k, 1}, sel22{k, 2}, w.outer_diameter * 1e3));
    L22{end + 1} = sprintf('%s %s: steady conductor %.1f degC, after 60 s %.1f degC', sel22{k, 1}, sel22{k, 2}, r.T(end, 1), interp1(r.t, r.T(:, 1), 60)); %#ok<SAGROW>
end
xlabel('Time [s]'); ylabel(ct_lbl(['Conductor temperature [' deg 'C]'])); ct_legend('southeast');
title(ct_lbl(sprintf('22 AWG wires, 8 A step, still air %g %sC', AMBIENT, deg)));
ct_save_fig(f, fullfile(out, 'fig25_22awg_insulation_transient.png'));
fid = fopen(fullfile(out, 'comparisons.txt'), 'a'); fprintf(fid, '%s\n', L22{:}); fclose(fid);
fprintf('%s\n', L22{:});
fprintf('results written to %s\n', out);
