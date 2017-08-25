function da1xdt = a1ODE(t,a1x)
  da1xdt = zeros(size(a1x));
  da1xdt(2) = a1x(3);
  da1xdt(3) = -(a1x(2)-3.5);
end