A = hilb(5);
[Q, R] = qr(A)
[AStar, gamma] = qr_decomposition(A);
Q2 = generate_q(AStar, gamma)
R2 = triu(AStar)