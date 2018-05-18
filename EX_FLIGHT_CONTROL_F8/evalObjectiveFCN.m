function J = evalObjectiveFCN(u,x,xref,Q,R,Ru)
% [J,Js,Ju]
%% Cost function of nonlinear MPC for dynamical system
%
% Inputs:
%   u:      time series of control input 
%   x:      time series of state
%   Ts:     controller sample time
%   xref:   time series of state references
%    Q:      State weights
%    R:      Penalization/weights on rate of change in u, uk - uk-1 
%    Ru:     Control input u penalization/weights
%
% Output:
%   J:      objective function cost
%

%% Cost Calculation
% Set initial cost and input
N = size(x,2);
J = zeros(N,1);
u0 = 0;

% Loop through each time step.
for ct=1:N
    
    % Accumulate state tracking cost
    J(ct) = abs(x(:,ct)-xref(:,ct))'*Q*abs(x(:,ct)-xref(:,ct));
    
    % Accumulate rate of change cost + cost of input
    if ct==1
        J(ct) = J(ct) + abs(u(:,ct)-u0)'*R*abs(u(:,ct)-u0) + abs(u(:,ct))'*Ru*abs(u(:,ct));
    else
        J(ct) = J(ct) + abs(u(:,ct)-u(:,ct-1))'*R*abs(u(:,ct)-u(:,ct-1)) + abs(u(:,ct))'*Ru*abs(u(:,ct));
    end
end
