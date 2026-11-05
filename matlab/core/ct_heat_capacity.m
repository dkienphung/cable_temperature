function C = ct_heat_capacity(m, T)
%CT_HEAT_CAPACITY  C_i = sum_m rho_m S_i,m c_p,m(T_i)  [J/(m K)] (left-hand sides of Eqs. 5-23).
C = zeros(1, m.n);
for k = 1:3
    C = C + m.rho(k) * m.S(:, k)' .* ct_mat_cp(m.mats{k}, T);
end
C = m.opt.capacity_scale * C;
end
