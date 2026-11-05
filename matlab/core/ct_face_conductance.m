function G = ct_face_conductance(m, T)
%CT_FACE_CONDUCTANCE  G_s = k_m(T_face) p_s / dx_m  [W/(m K)]  (the k p/dx terms of Eqs. 5-27).
if strcmp(m.opt.face_k, 'mean'), Tf = 0.5 * (T(1:end - 1) + T(2:end)); else, Tf = T(1:end - 1); end
k = zeros(1, m.n - 1);
for j = 1:3
    sel = m.seg_mat == j;
    k(sel) = ct_mat_k(m.mats{j}, Tf(sel));
end
G = k .* m.face_p ./ m.seg_dx;
end
