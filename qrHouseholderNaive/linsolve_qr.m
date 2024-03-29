%Solve a linear system Ax = b using naive QR decomposition

function x = linsolve_qr(A, b)
  [Q, R] = qr_decomposition(A);
  c = Q'*b;
  x = backsubs(R, c);
end