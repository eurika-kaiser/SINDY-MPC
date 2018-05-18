N  = 24;%24;%18;%12;%13;                        % Prediction horizon (number of iterations)
Nu  = N;                            % Control horizon (number of iterations)

% Q = [-1 0 0 1 5];                    % State weights
Q = [1 0 1 0 0];
R = 0;                              % du weights
Ru = 1;                             % u weights
LB = 0*ones(Nu,1);                  % Lower bound of control input
UB = 1*ones(Nu,1);                % Upper bound of control input
LBdu = nan;                         % Lower bound of control input
UBdu = nan;                         % Upper bound of control input
LBo = 0;                            % Lower bound of output
UBo = nan;                          % Upper bound of output
