% Generates Q from the matrix returned by qr_decomposition

function Q = generate_q(A, gamma)
  n = size(A)(1);
  m = size(A)(2);
  assert(length(gamma) == m, "Length of gamma and #columns of A differ");
  Q = eye(n);
  %Find the column index of the last reflector
  if n == m
    endindex = m-1;
  else
    endindex = m;
  endif
  %Apply reflectors Q1Q2 ... Qendindex, right to left
  for k = endindex:-1:1
    u = ones(n-k+1, 1);
    u(2:end, 1) = A(k+1:n, k);
    vT = gamma(k)*u';
    vT = vT*Q(k:end, k:end);
    Q(k:end, k:end) = Q(k:end, k:end) - u*vT;
  endfor
end      