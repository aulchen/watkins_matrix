%Given a real square matrix A, compute the QR decomposition

function [Q, A] = qr_decomposition_naive(A)
  n = size(A)(1);
  Q = eye(n);
  if n == 1
    return
  endif
  for k = 1:n-1
    %Compute the reflector
    [tau, gamma, u] = reflector(A(k:n, k));
    %Compute QkA
    vT = gamma*u';
    vT = vT*(A(k:n, k+1:n));
    A(k:n, k+1:n) = A(k:n, k+1:n) - u * vT;
    A(k, k) = -tau;
    A(k+1:n, k) = zeros(n-k, 1);
    %Compute Q1...Qk by right multiplying
    %Update the lower-right block
    v = gamma*u;
    v = Q(k:n, k:n)*v;
    Q(k:n, k:n) = Q(k:n, k:n)-(v*u');
    %Update the upper-right block
    if k > 1
      v = gamma*u;
      v = Q(1:k-1, k:n)*v;
      Q(1:k-1, k:n) = Q(1:k-1, k:n)-(v*u');
    endif
  endfor
end