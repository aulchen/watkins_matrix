% Returns the assembled Q for a given vector x

function Q = reflector_q(x)
  [tau, gamma, u] = reflector(x);
  Q = eye(length(x)) - u*(gamma*u');
end      