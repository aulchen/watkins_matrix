%Solve a linear system Ax = b using naive QR decomposition

function x = linsolve_qr_naive(A, b)
  [Q, R] = qr_decomposition_naive(A);
  c = Q'*b;
  x = backsubs(R, c);
end