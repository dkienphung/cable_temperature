function out = ct_iec60949(mode, a, b, S_mm2, theta_i, K, beta)
%CT_IEC60949  Adiabatic short-circuit temperature rise, I^2 t = K^2 S^2 ln((theta_f+beta)/(theta_i+beta)).
%   theta_f = ct_iec60949('final',   I, t, S_mm2, theta_i, K, beta)
%   I       = ct_iec60949('current', theta_f, t, S_mm2, theta_i, K, beta)
%   K       = ct_iec60949('K', sigma_c, rho20, [], [], [], beta)  [A s^0.5/mm2]
%   Copper: sigma_c = 3.45e6 J/(K m3), beta = 234.5 K, rho20 = 1.7241e-8 -> K = 226.
if nargin < 6 || isempty(K), K = 226; end
if nargin < 7 || isempty(beta), beta = 234.5; end
switch mode
    case 'final'
        out = (theta_i + beta) * exp(a.^2 .* b / (K^2 * S_mm2^2)) - beta;
    case 'current'
        out = K * S_mm2 * sqrt(log((a + beta) / (theta_i + beta)) / b);
    case 'K'
        out = sqrt(a * (beta + 20) / b) * 1e-6;
    otherwise
        error('ct_iec60949:mode', 'mode must be final, current or K');
end
end
