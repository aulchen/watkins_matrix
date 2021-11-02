A = [2 4 2 3;
     -2 -5 -3 -2;
     4 7 6 8;
     6 10 1 12]
b = [-3 3 -1 -16]';
x = linsolve(A, b)
x2 = linsolve_qr(A, b)