N  = 13;%13;15                        % Prediction horizon (number of iterations)
Nu  = N;                        % Control horizon (number of iterations)

Q = [25 0 0];                      % State weights
R = 0.05;%.1;                           % du weights
Ru = 0.05;%.1;                            % u weights
LB = -0.3*ones(Nu,1);        % Lower bound of control input, 0.3
UB = 0.5*ones(Nu,1);         % Upper bound of control input, 0.5
LBdu = nan;        % Lower bound of control input
UBdu = nan;         % Upper bound of control input
LBo = [-0.2,-1,-1];                 % Lower bound of output
UBo = [0.4,1,1];                  % Upper bound of output

% Reference trajectory
time_control = 0:dt:Duration;
r = 0.4*(-0.5./(1+exp(time_control./0.1-8)) + 1./(1+exp(time_control./0.1-30)) - 0.4);
figure,plot(time_control,r,'-k'), hold on

% xrefFUN = @(t) 0.4*(-0.5./(1+exp(t-8)) + 1./(1+exp(t-30)) - 0.4);

xrefFUN = @(t) 0.4*(-0.5./(1+exp(t/0.1-8)) + 1./(1+exp(t/0.1-30)) - 0.4);
% xrefFUN = @(t) 0.4*(-0.5./(1+exp(t-8/0.1)) + 1./(1+exp(t-30/0.1)) - 0.4);
plot(time_control,xrefFUN(time_control),'--r')



%% works
% N  = 13;%13;                        % Prediction horizon (number of iterations)
% Nu  = N;                        % Control horizon (number of iterations)
% 
% Q = [25 0 0];                      % State weights
% R = 0.;                           % du weights
% Ru = 0.1;                            % u weights
% LB = -0.3*ones(Nu,1);        % Lower bound of control input, 0.3
% UB = 0.5*ones(Nu,1);         % Upper bound of control input, 0.5
% LBdu = -0.1;        % Lower bound of control input
% UBdu = 0.1;         % Upper bound of control input
% LBo = [-0.2,-1,-1];                 % Lower bound of output
% UBo = [0.4,1,1];                  % Upper bound of output
% 
% % Reference trajectory
% time_control = 0:dt:Duration;
% r = 0.4*(-0.5./(1+exp(time_control./0.1-8)) + 1./(1+exp(time_control./0.1-30)) - 0.4);
% figure,plot(time_control,r,'-k'), hold on
% 
% % xrefFUN = @(t) 0.4*(-0.5./(1+exp(t-8)) + 1./(1+exp(t-30)) - 0.4);
% 
% xrefFUN = @(t) 0.4*(-0.5./(1+exp(t/0.1-8)) + 1./(1+exp(t/0.1-30)) - 0.4);
% % xrefFUN = @(t) 0.4*(-0.5./(1+exp(t-8/0.1)) + 1./(1+exp(t-30/0.1)) - 0.4);
% plot(time_control,xrefFUN(time_control),'--r')