A = [2 1 -1 3; -2 0 0 0; 4 1 -2 6; -6 -1 2 -3]
b = [13 -2 24 -14]';
x = linsolve(A, b)
x2 = linsolve_qr_naive(A, b)