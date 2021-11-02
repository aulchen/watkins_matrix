%Given a real vector x, computes the reflector that makes
%all entries except the first 0.

function [tau, gamma, u] = reflector(x)
  tau = sqrt(sum(x.^2));
  if x(1) < 0
    tau = -tau;
  endif
  u = x;
  u(1) = x(1) + tau;
  gamma = 2 / sum(u.^2);
end
