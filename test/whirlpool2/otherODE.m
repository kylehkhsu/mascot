function [t phi] = otherODE(x,u,tau)
  [t phi] = ode45(@bODE, [0 tau], x, [], u);
end



function dbxdt = bODE(t, bx,u)
  dbxdt = zeros(size(bx));
  dbxdt(1) = (-0.5 * bx(1) + 1 * bx(2)) * u(2);
  dbxdt(2) = (-0.5 * bx(1) + 0.5 * bx(2) + u(1)) * u(2);
%    dbxdt = dbxdt + d';
end