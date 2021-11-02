%Given a real square matrix A, compute the QR decomposition
%Outputs a matrix A containing R in the upper-right entries
%and u^k, except for the leading 1
%Also outputs a nx1 vector gamma

function [A, gamma] = qr_decomposition(A)
  n = size(A)(1);
  gamma = zeros(1, n);
  if n == 1
    return
  endif
  for k = 1:n-1
    %Compute the reflector
    [tau, gamma(k), u] = reflector(A(k:n, k));
    %Compute QkA
    vT = gamma(k)*u';
    vT = vT*(A(k:n, k+1:n));
    A(k:n, k+1:n) = A(k:n, k+1:n) - u * vT;
    %Update r
    A(k, k) = -tau;
    %Store u^k in the empty spaces of A
    A(k+1:n, k) = u(2:end);
  endfor
  gamma(n) = A(n, n);
end