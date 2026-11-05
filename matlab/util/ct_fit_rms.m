function v = ct_fit_rms(c, p, figkey)
%CT_FIT_RMS  RMS difference to the paper curves for heat-capacity / h multipliers p = [cap, h].
r = ct_simulate(ct_model(c, [], 'capacity_scale', abs(p(1)), 'h_scale', abs(p(2))), c.test_steps);
s = ct_paper_compare(r, figkey);
v = s.all_rms;
end
