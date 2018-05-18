function dx = HIVsys_ZURAKOWSKI(t,x,u,p)

try
    run_HIV_params
    
    dx = zeros(5,1);
    dx(1) = lambda1 - d*x(1) - alpha1*(1-eta*u(1))*x(1)*x(2);
    dx(2) = alpha1*(1-eta*u(1))*x(2)*x(1) - a*x(2) - p1*x(4)*x(2) - p2*x(5)*x(2);
    dx(3) = c2*x(1)*x(2)*x(3) - c2*q*x(2)*x(3) - b2*x(3);
    dx(4) = c1*x(2)*x(4) - b1*x(4);
    dx(5) = c2*q*x(2)*x(3) - h*x(5);
catch
    keyboard
end
