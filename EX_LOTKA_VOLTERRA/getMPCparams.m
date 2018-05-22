Q = [1 1];                      % State weights
R = 0.5;                        % du weights
Ru = 0.5;                       % u weights
B = [0; 1];
C = eye(Nvar);
D = 0;
LB = -20*ones(N,1);        % Lower bound of control input
UB = 20*ones(N,1);         % Upper bound of control input