function u = mynoise(t,tspan)
rng(0,'twister');
Nt = length(tspan);
uall = randn(1,Nt);
TF = t<=tspan;
u = uall(find(TF==1, 1, 'first'));
