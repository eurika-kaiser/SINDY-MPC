% MPC parameters
N  = 13;                        % Prediction horizon (number of iterations)
Nu  = N;                        % Control horizon (number of iterations)

Q = [25 0 0];                   % State weights
R = 0.05;                       % du weights, 0
Ru = 0.05;                      % u weights, 0.1
LB = -0.3*ones(Nu,1);           % Lower bound of control input, -0.3
UB = 0.5*ones(Nu,1);            % Upper bound of control input, 0.5
LBdu = nan;                     % Lower bound of control input rate, -0.1
UBdu = nan;                     % Upper bound of control input rate, 0.1
LBo = [-0.2,-1,-1];             % Lower bound of output
UBo = [0.4,1,1];                % Upper bound of output

% Reference trajectory
time_control = 0:dt:Duration;
r = 0.4*(-0.5./(1+exp(time_control./0.1-8)) + 1./(1+exp(time_control./0.1-30)) - 0.4);
figure,plot(time_control,r,'-k'), hold on

xrefFUN = @(t) 0.4*(-0.5./(1+exp(t/0.1-8)) + 1./(1+exp(t/0.1-30)) - 0.4);
plot(time_control,xrefFUN(time_control),'--r')