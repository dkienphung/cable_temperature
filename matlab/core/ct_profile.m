function [f, t_end] = ct_profile(spec)
%CT_PROFILE  Current / ambient profile as a function handle.
%   spec: scalar, function handle f(t), or step table [value duration; ...] (paper Tables 2 and 4).
t_end = [];
if isa(spec, 'function_handle')
    f = spec;
elseif isscalar(spec)
    f = @(t) spec + 0 * t;
elseif size(spec, 2) == 2
    edges = [0; cumsum(spec(:, 2))];
    vals = spec(:, 1);
    f = @(t) vals(min(max(sum(edges <= t), 1), numel(vals)));
    t_end = edges(end);
else
    error('ct_profile:spec', 'profile must be a scalar, a function handle or an [value duration] table');
end
end
