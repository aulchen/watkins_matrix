function X = qr_rect_test(x)
  A = [1 2; 3 4; 5 6];
  [Q, R] = qr(A)
  [Q2, R2] = qr_decomposition_rect(A)
end