%Given a real vector x, computes the reflector that makes
%all entries except the first 0.

function [tau, gamma, u] = reflector(x)
  n = length(x);
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
end
