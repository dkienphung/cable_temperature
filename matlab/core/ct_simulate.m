function r = ct_simulate(m, current, ambient, t_end, T_init, dt, record_every)
%CT_SIMULATE  Transient / dynamic simulation with the paper's implicit FDM + TDMA (Eqs. 3-27).
%   r = ct_simulate(m, current, ambient, t_end, T_init, dt, record_every)
%   current / ambient : scalar, function handle f(t) or [value duration] step table
%   ambient = [] -> cable.ambient_temperature;  t_end = [] -> length of the step table
%   T_init  = [] -> uniform ambient;  dt = [] -> m.opt.dt
%   Joule source, properties and radiation at time level j; convection implicit (Eq. 27);
%   the current of a step is evaluated at mid-step, the ambient at the end of the step.
%   Output r: t, T (nt x n), I, Ta, joule, p_conv, p_rad, stored [W/m], probes, x.
if nargin < 3 || isempty(ambient), ambient = m.cable.ambient_temperature; end
if nargin < 6 || isempty(dt), dt = m.opt.dt; end
if nargin < 7 || isempty(record_every), record_every = 1; end
[If, t_tab] = ct_profile(current);
Taf = ct_profile(ambient);
if nargin < 4 || isempty(t_end)
    if isempty(t_tab), error('ct_simulate:t_end', 't_end required'); end
    t_end = t_tab;
end
nsteps = round(t_end / dt);
if nargin < 5 || isempty(T_init), T = Taf(0) * ones(1, m.n);
elseif isscalar(T_init), T = T_init * ones(1, m.n);
else, T = T_init(:)';
end
nrec = floor(nsteps / record_every) + 2;
r.t = nan(nrec, 1); r.T = nan(nrec, m.n); r.I = nan(nrec, 1); r.Ta = nan(nrec, 1);
r.joule = nan(nrec, 1); r.p_conv = nan(nrec, 1); r.p_rad = nan(nrec, 1); r.stored = nan(nrec, 1);
k = 1; r.t(1) = 0; r.T(1, :) = T; r.I(1) = If(0); r.Ta(1) = Taf(0);
t = 0;
for j = 1:nsteps
    I = If(t + 0.5 * dt);
    Ta = Taf(t + dt);
    [a, b, cc, d, C, q, h, prad] = ct_coefficients(m, T, I, Ta, dt);
    Tn = ct_tdma(a, b, cc, d);
    if strcmp(m.opt.radiation, 'linearized')
        [~, s] = ct_surface_rad(m, T(end), Ta);
        prad = prad + s * (Tn(end) - T(end));
    end
    t = j * dt;
    if mod(j, record_every) == 0 || j == nsteps
        k = k + 1;
        r.t(k) = t; r.T(k, :) = Tn; r.I(k) = I; r.Ta(k) = Ta;
        r.joule(k) = sum(q); r.p_conv(k) = h * m.p_out * (Tn(end) - Ta); r.p_rad(k) = prad;
        r.stored(k) = sum(C .* (Tn - T)) / dt;
    end
    T = Tn;
end
f = {'t', 'T', 'I', 'Ta', 'joule', 'p_conv', 'p_rad', 'stored'};
for i = 1:numel(f), r.(f{i}) = r.(f{i})(1:k, :); end
r.probes = m.probes; r.x = m.x;
end
