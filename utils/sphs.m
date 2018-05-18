function u = sphs(Pf,K,t)
% Schroeder-phased harmonic sequence
% Pf = 1; % Fundamental period
% K = 10;
% u = zeros(size(time_test)); 
u = zeros(size(t));
for i = 1:K
    theta = 2*pi/K * sum([1:i]);
    u = u + sqrt(2/(K))*cos(2*pi * i*t/Pf + theta) ;
end
