%Given a real nxm, n>m matrix A, compute the QR decomposition

function [Q, A] = qr_decomposition_rect(A)
  n = size(A)(1); %Rows of A
  m = size(A)(2); %Columns of A
  Q = eye(n);
  if n == 1
    return
  endif
  for k = 1:m
    %Compute the reflector
    [tau, gamma, u] = reflector(A(k:n, k));
    %Compute QkA
    vT = gamma*u';
    vT = vT*(A(k:end, k+1:end));
    A(k:end, k+1:end) = A(k:end, k+1:end) - u * vT;
    A(k, k) = -tau;
    A(k+1:end, k) = zeros(n-k, 1);
    %Compute Q1...Qk by right multiplying
    %Update the lower-right block
    v = gamma*u;
    v = Q(k:end, k:end)*v;
    Q(k:end, k:end) = Q(k:end, k:end)-(v*u');
    %Update the upper-right block
    if k > 1
      v = gamma*u;
      v = Q(1:k-1, k:end)*v;
      Q(1:k-1, k:end) = Q(1:k-1, k:end)-(v*u');
    endif
  endfor
end