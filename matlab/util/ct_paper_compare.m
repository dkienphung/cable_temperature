function s = ct_paper_compare(r, figkey)
%CT_PAPER_COMPARE  Model vs digitized paper curves (Fig. 9 'fig9_cable70', Fig. 11 'fig11_cable150').
%   Returns RMS / mean / end values per measuring point (model - paper, K).
d = ct_json('paper_figures_digitized.json');
F = d.(figkey);
probes = {'conductor_insulation', 'insulation_jacket', 'jacket_air'};
keys = {'conductor_insulation', 'insulation_jacket', 'jacket_air_measured'};
allerr = [];
for i = 1:3
    P = F.(keys{i});
    e = interp1(r.t, ct_probe(r, probes{i}), P(:, 1)) - P(:, 2);
    s.(probes{i}) = struct('rms', sqrt(mean(e.^2)), 'mean', mean(e), 'max_abs', max(abs(e)), ...
        'end_model', r.T(end, r.probes.(probes{i})), 'end_paper', P(end, 2), 'paper', P);
    allerr = [allerr; e]; %#ok<AGROW>
end
s.all_rms = sqrt(mean(allerr.^2));
s.all_mean = mean(allerr);
end
