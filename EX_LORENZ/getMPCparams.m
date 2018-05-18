Q = [1 1 1];                      % State weights
R = 0.001;%0.5;                   % du weights
Ru = 0.001;%0.5;                  % u weights
% B = [1; 0; 0];
C = eye(Nvar);
D = 0;
LB = -50*ones(N,1);        % Lower bound of control input
UB = 50*ones(N,1);         % Upper bound of control input