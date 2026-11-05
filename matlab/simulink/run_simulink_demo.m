%RUN_SIMULINK_DEMO  Build the Simulink model, simulate the 70 mm2 test sequence of the paper and
%   compare it with ct_simulate (both use the same discrete equations).  Requires Simulink.
root = setup_paths();
cable = ct_cable_preset('cable70');
mdl = build_cable_simulink_model(cable);
out = sim(mdl);
Tp = out.get('T_probes');                  % [centre, cond-ins, ins-jacket, jacket-air] at t = 0, dt, ...
r = ct_simulate(ct_model(cable), cable.test_steps);
k = min(size(Tp, 1), numel(r.t));
fprintf('max |T_simulink - T_matlab| = %.2e K\n', max(max(abs(Tp(1:k, :) - r.T(1:k, [1 10 19 28])))));
figure; plot(r.t(1:k), Tp(1:k, 2:4)); grid on;
xlabel('Time [s]'); ylabel('Temperature [degC]'); legend('conductor-insulation', 'insulation-jacket', 'jacket-air');
title('Simulink model of the 70 mm2 cable (paper Table 2)');
