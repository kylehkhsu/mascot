function phi = RK(x, u, nint, tau)
  h = tau/nint;
  dim = size(x,2);
  k = zeros(4,dim);
  tmp = zeros(size(x));
  for t = 1:nint
    k(1,:) = bODE(x, u);
    for i = 1:dim
      tmp(i) = x(i) + h/2 * k(1,i);    
    end
    k(2,:) = bODE(tmp, u);
    for i = 1:dim
      tmp(i) = x(i) + h/2 * k(2,i);
    end
    k(3,:) = bODE(tmp, u);
    for i = 1:dim
      tmp(i) = x(i) + h * k(3,i);
    end
    k(4,:) = bODE(tmp, u);
    for i = 1:dim
      x(i) = x(i) + h/6 * (k(1,i) + 2*k(2,i) + 2*k(3,i) + k(4,i));
    end
  end
  phi = x;
end

function dbxdt = bODE(bx,u)
  dbxdt = zeros(size(bx));
  dbxdt(1) = (-0.5 * bx(1) + 1 * bx(2)) * u(2);
  dbxdt(2) = (-0.5 * bx(1) + 0.5 * bx(2) + u(1)) * u(2);
%    dbxdt = dbxdt + d';
end