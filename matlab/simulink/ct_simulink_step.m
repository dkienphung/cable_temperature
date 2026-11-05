function Tn = ct_simulink_step(T, I, Ta)
%CT_SIMULINK_STEP  One time step of the paper's FDM/TDMA model, called (extrinsically) from the
%   MATLAB Function block of the Simulink model built by build_cable_simulink_model.
%   The model struct is read once from the base-workspace variable CT_MODEL (set by the builder);
%   run "clear ct_simulink_step" after changing CT_MODEL.
persistent m
if isempty(m), m = evalin('base', 'CT_MODEL'); end
[a, b, c, d] = ct_coefficients(m, T(:)', I, Ta, m.opt.dt);
Tn = ct_tdma(a, b, c, d)';
end
