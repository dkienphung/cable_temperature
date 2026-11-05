function h = ct_h_surface(Ts, Ta, D, H, model, V, angle_deg)
%CT_H_SURFACE  Convective coefficient at the jacket-air node: natural convection, or (V > 0)
%   the larger of natural and IEEE 738 forced convection (as IEEE 738 prescribes).
if nargin < 6 || isempty(V), V = 0; end
if nargin < 7 || isempty(angle_deg), angle_deg = 90; end
h = ct_h_natural(Ts, Ta, D, H, model);
if V > 0
    h = max(h, ct_h_forced(Ts, Ta, D, V, H, angle_deg));
end
end
