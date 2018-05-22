function J = evalObjectiveFCN(u,x,xref,Q,R,Ru)
% [J,Js,Ju]
%% Cost function of nonlinear MPC for dynamical system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   xref:   state references, constant from time k+1 to k+N
%
% Output:
%   J:      objective function cost
%

%% Nonlinear MPC design parameters

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
        J(ct) = J(ct) + abs(u(:,ct)-u0)'*R*abs(u(:,ct)-u0) + abs(u(:,ct))'*Ru*abs(u(:,ct));
    else
        J(ct) = J(ct) + abs(u(:,ct)-u(:,ct-1))'*R*abs(u(:,ct)-u(:,ct-1)) + abs(u(:,ct))'*Ru*abs(u(:,ct));
    end
end
