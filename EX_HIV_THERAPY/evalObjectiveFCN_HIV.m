function J = evalObjectiveFCN_HIV(u,x,xref,Q,R,Ru)

%% Cost function of nonlinear MPC for HIV model
%
% Inputs:
%   u:      time series of control input 
%   x:      time series of state
%   xref:   time series of state references
%    Q:      State weights
%    R:      Penalization/weights on rate of change in u, uk - uk-1 
%    Ru:     Control input u penalization/weights
%
% Output:
%   J:      objective function cost
%

%% Cost Calculation
% Set initial cost and input.
N = size(x,2);
J = zeros(N,1);
u0 = 0;

% Loop through each time step.
for ct=1:N
    
    % Accumulate state tracking cost from x(k+1) to x(k+N).
    J(ct) = abs(x(:,ct)-xref(:,ct))'*Q*abs(x(:,ct)-xref(:,ct));
    
    % accumulate MV rate of change cost from u(k) to u(k+N-1).
    if ct==1
        J(ct) = J(ct) + abs(u(:,ct)-u0)'*R*abs(u(:,ct)-u0) + Ru*abs(ct);
    else
        J(ct) = J(ct) + abs(u(:,ct)-u(:,ct-1))'*R*abs(u(:,ct)-u(:,ct-1)) + Ru*abs(ct);
    end
end
