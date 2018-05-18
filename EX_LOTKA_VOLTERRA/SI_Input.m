function u = SI_Input(t,A,Pf,K,tspan,taulim, Nswitch, states)

if t<=20
    u = A(1)*chirp(t,[],20,1.).^2;
elseif 20<t && t<=50
    u = A(2)*sphs(Pf,K,t);
elseif 50<t && t<=60
    u = A(3)*mynoise(t,tspan);
elseif 60<t && t<=100
    Toffset = 60;
    u = A(4)*prbs(taulim, Nswitch, states, t,Toffset);
end