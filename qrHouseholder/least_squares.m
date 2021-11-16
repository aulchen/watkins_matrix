%Given overdetermined linear system Ax = b with A being nxm and having full rank,
%compute the least squares solution using Householder reflectors.
%Q is nxm, R is mxm, and b is nx1.

function [x, residual] = least_squares(A, b)
  n = size(A)(1); %Number of rows
  m = size(A)(2); %Number of columns
  [A, gamma] = qr_decomposition(A);
  R = triu(A);
  %Repeatedly apply reflectors Q1, Q2... Qn-1 to b
  %Since QTb = Qm-1...Q1b = c
  for k = 1:m %For every column
    %Extract u from A
    u = ones(n-k+1, 1);
    u(2:end, 1) = A(k+1:end, k);
    %Compute Qkb using reflectors
    vT = gamma(k)*u';
    vT = vT*b(k:end);
    b(k:end) = b(k:end) - u*vT;
  endfor
  %Compute the residual as the norm of the m+1 to nth terms of b
  residual = norm(b(m+1:end), 2);
  %Backsolve Rx = c_hat, where c_hat is the 1 to mth terms of b
  x = backsubs(R, b(1: m));
end