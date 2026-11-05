function o = ct_options(varargin)
%CT_OPTIONS  Model options (defaults = the paper's model, corrected Eq. 25).
%   nodes          [NA NB NC]        nodes per material (paper: 9 9 9)
%   dt             10 s              time step (paper)
%   convection     'ieee738'         'ieee738' (Eq. 25 corrected) | 'paper' (as printed) | 'churchill_chu'
%   elevation      []                H of Eq. 26 [m]; [] -> cable.elevation
%   air_speed      0                 indoor air speed [m/s]; > 0 adds IEEE 738 forced convection
%   air_angle      90                angle air flow / cable axis [deg]
%   radiation      'explicit'        'explicit' (paper Eq. 23) | 'linearized'
%   joule_area     'effective'       'effective' (r = rho/S_eff) | 'geometric' (rho/(pi R_A^2))
%   face_k         'mean'            k on a face at the mean of the two node temperatures
%   convection_off, radiation_off, adiabatic : verification switches
%   joule_scale, h_scale, rad_scale, capacity_scale : sensitivity multipliers (1 = paper)
o = struct('nodes', [9 9 9], 'dt', 10, 'convection', 'ieee738', 'elevation', [], 'air_speed', 0, ...
    'air_angle', 90, 'radiation', 'explicit', 'joule_area', 'effective', 'face_k', 'mean', ...
    'convection_off', false, 'radiation_off', false, 'adiabatic', false, ...
    'joule_scale', 1, 'h_scale', 1, 'rad_scale', 1, 'capacity_scale', 1);
for i = 1:2:numel(varargin)
    if ~isfield(o, varargin{i}), error('ct_options:field', 'unknown option %s', varargin{i}); end
    o.(varargin{i}) = varargin{i + 1};
end
end
