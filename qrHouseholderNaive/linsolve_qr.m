%Solve a linear system Ax = b using QR decomposition
%b is assumed to be a column vector

function x = linsolve_qr(A, b)
  n = length(b);
  %Make b a column vector
  [A, gamma] = qr_decomposition(A);
  R = triu(A);
  %Repeatedly apply reflectors Q1, Q2... Qn-1 to b
  %Since QTb = Qn-1...Q1b = c
  for k = 1:n-1
    %Extract u from A
    u = ones(n-k+1, 1);
    u(2:end, 1) = A(k+1:n, k);
    %Compute Qkb using reflectors
    vT = gamma(k)*u';
    vT = vT*b(k:end);
    b(k:end) = b(k:end) - u*vT;
  endfor
  %Backsolve Rx = Qb
  x = backsubs(A, b);
end