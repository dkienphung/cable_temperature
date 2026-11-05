function mdl = build_cable_simulink_model(cable, mdl, opt)
%BUILD_CABLE_SIMULINK_MODEL  Build a Simulink model of the paper's cable thermal model.
%
%   build_cable_simulink_model                        % cable70, Table 2 current steps
%   build_cable_simulink_model(ct_cable_preset('cable150'))
%   build_cable_simulink_model(cable, 'my_model', ct_options('convection', 'churchill_chu'))
%
%   Structure (fixed-step discrete, sample time dt = opt.dt):
%
%     Clock -> +dt/2 -> 1-D Lookup (current steps, flat) --I-->+-----------------+
%     Constant T_air ----------------------------------Ta-->  | MATLAB Function |--T(k+1)--> Unit Delay --T(k)--+
%                                                    +----->  | FDM/TDMA step   |                    (IC = T0)  |
%                                                    |        +-----------------+                               |
%                                                    +--------------------------------------------------------- +
%     T(k) -> Selector [centre, conductor-insulation, insulation-jacket, jacket-air] -> Scope, To Workspace
%
%   The MATLAB Function block calls ct_simulink_step (coder.extrinsic), i.e. exactly the same
%   coefficients (Eqs. 5-27) and TDMA as ct_simulate, with the current evaluated at mid-step
%   like ct_simulate.  Requires Simulink.  NOTE: written for MATLAB/Simulink R2019b+ and not
%   executed in the environment where it was written (no Simulink there) - see README.
if nargin < 1 || isempty(cable), cable = ct_cable_preset('cable70'); end
if nargin < 2 || isempty(mdl), mdl = 'cable_thermal_sim'; end
if nargin < 3 || isempty(opt), opt = ct_options(); end
m = ct_model(cable, opt);
n = m.n;
steps = cable.test_steps;
if isempty(steps), steps = [100 3600]; end
edges = [0; cumsum(steps(:, 2))];
% workspace variables used by the block parameters
assignin('base', 'CT_MODEL', m);
assignin('base', 'CT_DT', opt.dt);
assignin('base', 'CT_T0', cable.ambient_temperature * ones(n, 1));
assignin('base', 'CT_I_BP', edges(1:end - 1)');
assignin('base', 'CT_I_TAB', steps(:, 1)');
assignin('base', 'CT_TAIR', cable.ambient_temperature);
assignin('base', 'CT_PROBES', [m.probes.center, m.probes.conductor_insulation, m.probes.insulation_jacket, m.probes.jacket_air]);
clear ct_simulink_step

if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', 'CT_DT', 'StopTime', num2str(edges(end)));

add_block('simulink/Sources/Clock', [mdl '/Clock'], 'Position', [30 40 60 70]);
add_block('simulink/Math Operations/Bias', [mdl '/mid-step'], 'Bias', 'CT_DT/2', 'Position', [100 40 140 70]);
add_block('simulink/Lookup Tables/1-D Lookup Table', [mdl '/Current steps'], 'Position', [180 30 260 80]);
set_param([mdl '/Current steps'], 'Table', 'CT_I_TAB', 'BreakpointsForDimension1', 'CT_I_BP', ...
          'InterpMethod', 'Flat', 'ExtrapMethod', 'Clip');
add_block('simulink/Sources/Constant', [mdl '/T_air'], 'Value', 'CT_TAIR', 'Position', [180 120 260 150]);
add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/FDM step'], 'Position', [330 40 470 160]);
add_block('simulink/Discrete/Unit Delay', [mdl '/T_k'], 'InitialCondition', 'CT_T0', 'SampleTime', 'CT_DT', ...
          'Position', [530 85 570 115]);
add_block('simulink/Signal Routing/Selector', [mdl '/probes'], 'InputPortWidth', num2str(n), ...
          'Indices', 'CT_PROBES', 'Position', [620 85 660 115]);
add_block('simulink/Sinks/Scope', [mdl '/Temperatures'], 'Position', [720 70 760 130]);
add_block('simulink/Sinks/To Workspace', [mdl '/To Workspace'], 'VariableName', 'T_probes', ...
          'SaveFormat', 'Array', 'Position', [720 160 800 190]);
add_block('simulink/Sinks/To Workspace', [mdl '/All nodes'], 'VariableName', 'T_nodes', ...
          'SaveFormat', 'Array', 'Position', [720 220 800 250]);

% MATLAB Function block code
code = sprintf(['function Tn = fdm_step(T, I, Ta)\n' ...
                '%%#codegen\n' ...
                '%% One step of the paper''s FDM/TDMA model (see ct_simulink_step / ct_coefficients)\n' ...
                'coder.extrinsic(''ct_simulink_step'');\n' ...
                'Tn = zeros(%d, 1);\n' ...
                'Tn = ct_simulink_step(T, I, Ta);\n'], n);
rt = sfroot;
blk = rt.find('-isa', 'Stateflow.EMChart', 'Path', [mdl '/FDM step']);
blk.Script = code;

add_line(mdl, 'Clock/1', 'mid-step/1');
add_line(mdl, 'mid-step/1', 'Current steps/1');
add_line(mdl, 'T_k/1', 'FDM step/1', 'autorouting', 'on');
add_line(mdl, 'Current steps/1', 'FDM step/2');
add_line(mdl, 'T_air/1', 'FDM step/3');
add_line(mdl, 'FDM step/1', 'T_k/1');
add_line(mdl, 'T_k/1', 'probes/1', 'autorouting', 'on');
add_line(mdl, 'probes/1', 'Temperatures/1');
add_line(mdl, 'probes/1', 'To Workspace/1', 'autorouting', 'on');
add_line(mdl, 'T_k/1', 'All nodes/1', 'autorouting', 'on');
save_system(mdl);
fprintf('Simulink model %s.slx built (%d nodes, dt = %g s, stop time %g s)\n', mdl, n, opt.dt, edges(end));
end
