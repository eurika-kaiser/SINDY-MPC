function J = evalObjectiveFCN(u,x,xref,Q,R,Ru)
% [J,Js,Ju]
%% Cost function of nonlinear MPC for dynamical system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%   xref:   state references, constant from time k+1 to k+N
%   u0:     previous controller output at time k-1
%
% Output:
%   J:      objective function cost
%

%% Nonlinear MPC design parameters
% Q = diag([1,1,1]);
% R = 0.01;

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

%% Cost Calculation
% % Set initial cost and input.
% N = size(x,2);
% Js = zeros(N,1);
% Ju = zeros(N,1);
% u0 = 0;
% 
% % Loop through each time step.
% for ct=1:N
%     
%     % Accumulate state tracking cost from x(k+1) to x(k+N).
%     Js(ct) = (x(:,ct)-xref(:,ct))'*Q*(x(:,ct)-xref(:,ct));
%     
%     % accumulate MV rate of change cost from u(k) to u(k+N-1).
%     if ct==1
%         Ju(ct) = (u(:,ct)-u0)'*R*(u(:,ct)-u0) + u(:,ct)'*Ru*u(:,ct);
%     else
%         Ju(ct) = (u(:,ct)-u(:,ct-1))'*R*(u(:,ct)-u(:,ct-1)) + u(:,ct)'*Ru*u(:,ct);
%     end
% end