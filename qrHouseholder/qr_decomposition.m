%Given a full rank real square nxm, n >= m matrix A, compute the QR decomposition
%Outputs a matrix A containing R in the upper-right entries
%and u^(k) in the lower-right, except for the leading 1.
%Also outputs a mx1 vector gamma

function [A, gamma] = qr_decomposition(A)
  n = size(A)(1); %Number of rows
  m = size(A)(2); %Number of columns
  assert(n >= m, "Matrix does not have more rows than columns");
  gamma = zeros(1, m);
  if n == 1
    return
  endif
  if n == m
    endindex = n - 1; %If matrix is square, we need to avoid using last column
  else
    endindex = m;
  endif
  for k = 1:endindex
    %Compute the reflector
    [tau, gamma(k), u] = reflector(A(k:end, k));
    %If any gamma_k == 0, then A is not full rank.
    if gamma(k) == 0
      error("A is not full rank.");
    endif
    %Compute QkA
    vT = gamma(k)*u';
    vT = vT*(A(k:end, k+1:end));
    A(k:end, k+1:end) = A(k:end, k+1:end) - u * vT;
    %Update r
    A(k, k) = -tau;
    %Store u^k in the empty spaces of A
    A(k+1:end, k) = u(2:end);
  endfor
  if n == m %Do this if A is square
    if A(n, n) == 0
      error("A is not full rank.");
    endif
    gamma(n) = A(n, n);
  endif
end