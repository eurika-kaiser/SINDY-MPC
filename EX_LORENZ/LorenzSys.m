function xdot = LorenzSys(t,x,u,p) 


xdot = [p.SIGMA*(-x(1)+x(2))+u;
        p.RHO*x(1)-x(2)-x(1)*x(3);
        -p.BETA*x(3)+x(1)*x(2)];