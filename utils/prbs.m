function u = prbs(taulim, Nswitch, states, t,Toffset,seed)
if nargin==6
    rng(seed,'twister')
    
else
    rng(0,'twister');
end

tau = (taulim(2)-taulim(1))*rand(Nswitch,1) + taulim(1);
tau = cumsum(tau)+Toffset;
ub = states(randi(length(states),Nswitch,1));

try
    TF = t<tau;
    u = ub(find(TF==1, 1, 'first'));
catch
    keyboard
end

end