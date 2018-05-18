function dx = HIVsys_KWON(t,x,u,p)

try
    run_HIV_params
    
    dx = zeros(5,1);
    dx(1) = lambda1 - d*x(1) - (1-u(1))*k1*x(4)*x(1) - k2*x(5)*x(1);
    dx(2) = (1-u(1))*k1*(1-mu)*x(4)*x(1) - delta*x(2);
    dx(3) = (1-u(1))*k1*mu*x(4)*x(1) + k2*x(5)*x(1) - delta*x(3);
    dx(4) = alpha1*x(2) - c*x(4);
    dx(5) = alpha2*x(3) - c*x(5);
%     dx(1) = lambda1 - d*x(1) - k1*x(1)*x(4) - k2*x(1)*x(5) + k1*x(1)*x(4)*u(1);
%     dx(2) =  - delta*x(2) + k1*(1-mu)*x(1)*x(4) - k1*(1-mu)*x(1)*x(4)*u(1);
%     dx(3) =  - delta*x(3) + k1*mu*x(1)*x(4) + k2*x(1)*x(5) - k1*mu*x(1)*x(4)*u(1);
%     dx(4) = alpha1*x(2) - c*x(4);
%     dx(5) = alpha2*x(3) - c*x(5);
catch
    keyboard
end
