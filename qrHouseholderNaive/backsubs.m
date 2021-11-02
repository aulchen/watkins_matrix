%Solve Ax = b using backwards substitution, where A is upper triangular
function x = backsubs(A, b)
  n = length(b);
  x = b;
  for i = n:-1:1
    for j = i+1:n
      x(i) = x(i) - A(i,j)*x(j);
    endfor
    assert(A(i, i) != 0, "Diagonal entry of U is 0");
    x(i) = x(i)/A(i, i);
  endfor
end