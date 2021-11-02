%Solve a linear system Ax = b using the QR decomposition

function x = linsolve_qr(A, b)
  [Q, R] = qr_decomposition_naive(A);
  c = Q'*b;
  x = backsubs(R, c);
end