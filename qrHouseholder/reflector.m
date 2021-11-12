%Given a real vector x, computes the reflector that makes
%all entries except the first 0.

function [tau, gamma, u] = reflector(x)
  n = length(x);
  beta = max(abs(x));
  if beta == 0
    tau = 0;
    gamma = 0;
    u = zeros(n, 1);
    return
  endif
  x = x / beta;
  tau = sqrt(sum(x.^2));
  if x(1) < 0
    tau = -tau;
  endif
  u = x;
  u(1) = x(1) + tau;
  %Rescale u
  u(2:n) = u(2:n) / u(1);
  gamma = u(1)/tau;
  u(1) = 1;
  tau = tau * beta;
end
