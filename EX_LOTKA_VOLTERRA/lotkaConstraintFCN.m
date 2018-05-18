function [c, ceq] = lotkaConstraintFCN(u,x,Ts,N,p)
%% Constraint function of nonlinear MPC for Lotka-Volterra system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%
% Output:
%   c:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints (empty)
%

%% Nonlinear MPC design parameters
% Predator population size always positive: >0, min pop size
zMin = [10];

%% Inequality constraints calculation
c = zeros(N,1);
% Apply N population size constraints across prediction horizon, from time
% k+1 to k+N
xk = x;
uk = u(1);
for ct=1:N
    % obtain new cart position at next prediction step
    xk1 = rk4u(@sparseGalerkinControl_Discrete,xk,uk,Ts,1,[],p); 
    
    % -z + zMin < 0 % lower bound
    c(ct) = -xk1(2)+zMin; %c(2*ct-1)
    % z - zMax < 0 % upper bound
    %c(2*ct) = xk1(1)-zMax;
    % update plant state and input for next step
    xk = xk1;
    if ct<N
        uk = u(ct+1);
    end
end
%% No equality constraints
ceq = [];

