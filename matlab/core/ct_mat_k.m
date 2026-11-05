function v = ct_mat_k(m, T)
%CT_MAT_K  Thermal conductivity k(T) [W/(m K)].
v = ct_interp_clamp(m.kT, m.kV, T);
end
