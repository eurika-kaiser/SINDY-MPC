function dy = lorenzcontrol_discrete(t,y,u,p)

dy = [p.SIGMA*(-y(1)+y(2))+u;
    p.RHO*y(1)-y(2)-y(1)*y(3);
    -p.BETA*y(3)+y(1)*y(2)];


