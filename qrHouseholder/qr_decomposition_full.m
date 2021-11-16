%Computes the QR decomposition of a nxm, n >= m full rank matrix A
%Outputs the matrices Q and R, where Q is orthogonal, R is upper triangular, and QR = A

function [Q, R] = qr_decomposition_full(A)
  [A, gamma] = qr_decomposition(A);
  R = triu(A);
  Q = generate_q(A, gamma);
end