function h = ct_h_natural(Ts, Ta, D, H, model)
%CT_H_NATURAL  Natural-convection coefficient h [W/(m2 K)] of a horizontal cylinder, diameter D [m].
%   model = 'paper'         : Eq. (25) exactly as printed, h = 3.645/pi rho^0.25 D^-0.25 dT^0.25
%   model = 'ieee738'       : Eq. (25) corrected (IEEE 738-2012 Eq. 5c), rho^0.5
%   model = 'churchill_chu' : Churchill & Chu (1975), valid 1e-5 <= Ra <= 1e12 (thin wires)
if nargin < 4 || isempty(H), H = 0; end
if nargin < 5, model = 'ieee738'; end
switch model
    case 'churchill_chu'
        h = ct_h_churchill_chu(Ts, Ta, D, H);
        return
    case 'paper'
        e = 0.25;
    case 'ieee738'
        e = 0.5;
    otherwise
        error('ct_h_natural:model', 'unknown convection model %s', model);
end
dT = abs(Ts - Ta);
rho = ct_air_density(Ts, Ta, H);
h = 3.645 / pi * rho.^e .* D.^(-0.25) .* dT.^0.25;
end
