function J = ObjectiveFCN(u,x,Ts,N,Nu,xref,u0,p,Q,R,Ru)
%% Cost function of nonlinear MPC for Lotka-Volterra system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%   Nu:     control horizon [not implemented]
%   xref:   state references, varying from time k+1 to k+N
%   u0:     previous controller output at time k-1
%    p:      Parameters for model
%    Q:      State weights
%    R:      Penalization/weights on rate of change in u, uk - uk-1 
%    Ru:     Control input u penalization/weights
%
% Output:
%   J:      objective function cost
%

%% Nonlinear MPC design parameters

%% Cost Calculation
% Set initial plant states, controller output and cost
xk = x;
uk = u(1);
J = 0;

% Loop through each prediction step
for ct=1:N

    % Obtain plant state at next prediction step
    xk1 = rk4u(@F8Sys,xk,uk,Ts,1,[],p);
    
    % Accumulate state tracking cost from x(k+1) to x(k+N)
    J = J + (xk1-xref(:,ct))'*Q*(xk1-xref(:,ct));
    
    % Accumulate rate of change cost from u(k) to u(k+N-1)
    if ct==1
        J = J + (uk-u0)'*R*(uk-u0) + uk'*Ru*uk;
    else
        J = J + (uk-u(ct-1))'*R*(uk-u(ct-1)) + uk'*Ru*uk;
    end
    
    % Update xk and uk for the next prediction step
    xk = xk1;
    if ct<N
        uk = u(ct+1);
    end
end


