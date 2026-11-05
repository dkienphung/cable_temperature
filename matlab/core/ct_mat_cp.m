function v = ct_mat_cp(m, T)
%CT_MAT_CP  Specific heat c_p(T) [J/(kg K)].
v = ct_interp_clamp(m.cpT, m.cpV, T);
end
