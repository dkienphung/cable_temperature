function T = ct_steady_state(m, I, Ta, T_guess)
%CT_STEADY_STATE  Steady state (dt -> Inf) of the same discrete equations, Picard iteration with
%   Newton-linearized radiation.  Returns the node temperatures [degC] (row vector).
if nargin < 3 || isempty(Ta), Ta = m.cable.ambient_temperature; end
if nargin < 4 || isempty(T_guess), T = (Ta + 10) * ones(1, m.n); else, T = T_guess(:)'; end
for it = 1:500
    [a, b, c, d] = ct_coefficients(m, T, I, Ta, Inf);
    Tn = ct_solve_tridiag(a, b, c, d);
    if max(abs(Tn - T)) < 1e-9 && all(isfinite(Tn)), T = Tn; return; end
    if it > 51, T = 0.5 * T + 0.5 * Tn; else, T = Tn; end
end
error('ct:steady', 'steady state did not converge');
end
