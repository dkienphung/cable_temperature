function q = ct_q_radiation(Ts, Ta, D, emissivity)
%CT_Q_RADIATION  Radiated power per unit length [W/m], paper Eq. (23): eps sigma pi D (Ts^4 - Ta^4).
if nargin < 4, emissivity = 0.85; end
c = ct_const();
q = emissivity * c.SIGMA * pi * D .* ((Ts + c.KELVIN).^4 - (Ta + c.KELVIN).^4);
end
