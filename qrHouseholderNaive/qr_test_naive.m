function X = qr_test(x)
  A = pascal(x);
  [Q, R] = qr(A)
  [Q2, R2] = qr_decomposition_naive(A)
end