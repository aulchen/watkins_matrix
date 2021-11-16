%Fit a quadratic polynomial for the data
%The data should fit phi(t) = 1 + t + t^2 exactly

t = [-1 -.75 -.5 0 .25 .5 .75]';
y = [1 .8125 .75 1 1.3125 1.75 2.3125]';
A = zeros(size(t)(1), 3); %length of t rows, 3 columns
A(:, 1) = ones(7, 1);
A(:, 2) = t;
A(:, 3) = t.^2;
[x, residual] = least_squares(A, y)