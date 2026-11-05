function x = ct_solve_tridiag(a, b, c, d)
%CT_SOLVE_TRIDIAG  Same tri-diagonal system as ct_tdma, solved with a sparse direct solver
%   (identical result to rounding; faster in Octave).  Used for the steady state.
n = numel(d);
A = sparse([1:n, 2:n, 1:n - 1], [1:n, 1:n - 1, 2:n], [b(:)', a(2:n), c(1:n - 1)], n, n);
x = (A \ d(:))';
end
