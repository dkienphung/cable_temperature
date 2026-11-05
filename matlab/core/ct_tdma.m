function x = ct_tdma(a, b, c, d)
%CT_TDMA  Thomas algorithm: a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i  (a_1 = c_n = 0), Eqs. (3)-(4).
n = numel(d);
cp = zeros(1, n); dp = zeros(1, n);
cp(1) = c(1) / b(1); dp(1) = d(1) / b(1);
for i = 2:n
    mm = b(i) - a(i) * cp(i - 1);
    if i < n, cp(i) = c(i) / mm; end
    dp(i) = (d(i) - a(i) * dp(i - 1)) / mm;
end
x = zeros(1, n);
x(n) = dp(n);
for i = n - 1:-1:1
    x(i) = dp(i) - cp(i) * x(i + 1);
end
end
