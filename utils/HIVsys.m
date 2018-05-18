function dx = HIVsys(t,x,u)

try
    lambda1 = 10; % [cells/(mm^3 day)]
    d1 = 0.01; %[1/day]
    k1 = 8e-4; %[mm^3/(virions day)]
    m1 = 0.01; %[mm^3/(cells day)]
    rho1 = 1; %[virions/cells]
    delta = 0.7; %[1/day]
    f = 0.34; %none
    lambdaE = 1e-3; %[cells/(mm^3 day)]
    bE = 0.3; %[1/day]
    Kb = 0.1; %[cells/mm^3]
    lambda2 = 31.98e-3; %[cells/(mm^3 day)]
    d2 = 0.01; %[1/day]
    k2 = 0.1; %[mm^3/(virions day)]
    m2 = 0.01; %[mm^3/(cells day)]
    rho2 = 1; %[virions/cells]
    c = 13; %[1/day]
    NT = 100; %[virions/day]
    deltaE = 0.1; %[1/day]
    dE = 0.25; %[1/day]
    Kd = 0.5; %[cells/mm^3]

    dx = zeros(6,1);
    dx(1) = lambda1 - d1*x(1) - (1-u(1))*k1*x(5)*x(1);
    dx(2) = lambda2 - d2*x(2) - (1-f*u(1))*k2*x(5)*x(2);
    dx(3) = (1-u(1))*k1*x(5)*x(1) - delta*x(3) - m1*x(6)*x(3);
    dx(4) = (1-f*u(1))*k2*x(5)*x(2) - delta*x(4) - m2*x(6)*x(4);
    dx(5) = (1-u(2))*NT*delta*(x(3)+x(4)) - (c+(1-u(1))*rho1*k1*x(1) + (1-f*u(1))*rho2*k2*x(2))*x(5);
    dx(6) = lambdaE + bE*( (x(3)+x(4))/(x(3)+x(4)+Kb) )*x(6) - dE*( (x(3)+x(4))/(x(3)+x(4)+Kd) )*x(6) - deltaE*x(6);
catch
    keyboard
end
% HIV model in Stengel, ...
% try
%     a1 = 2.4;
%     a2 = 2.4*10^(-5);
%     a3 = 1200;
%     a4 = 0.24;
%     a5 = 10;
%     a6 = 0.02;
%     a7 = 0.03;
%     a8 = 1500;
%     a9 = 0.003;
%     
%     u2 = 0;
%     u3 = 0;
%     u4 = 0;
%     
%     dx = zeros(4,1);
%     dx(1) = -a1*x(1)-a2*x(1)*x(2)*(1-u2)+a3*a4*x(4)*(1-u1);
%     dx(2) = a5/(1+x(1))-a2*x(1)*x(2)*(1-u2)*(1-u4)-a6*x(2)+a7*(1-(x(2)+x(3)+x(4))/a8)*x(2)*(1+u3);
%     dx(3) = a2*x(1)*x(2)*(1-u2)*(1-u4)-a9*x(3)-a6*x(3);
%     dx(4) = a9*x(3)-a4*x(4);
% catch
%     keyboard
% end


