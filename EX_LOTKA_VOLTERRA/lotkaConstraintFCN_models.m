function [c, ceq] = lotkaConstraintFCN_models(u,x,N,p,select_model)
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

%% Integrate system
if strcmp(select_model,'DelayDMDc')
     [xk,~] = lsim(p.sys,[p.udelay(1:N);u'],[0:N-1].*p.dt,[p.xdelay(:,1); x]-[p.xmean; p.xmean]);
    xk = xk(:,3:4);
    xk = xk + repmat(p.xmean',[N 1]); xk = xk';
elseif strcmp(select_model,'DMDc')
    [xk,~] = lsim(p.sys,[u',0],[0:N].*p.dt,x-p.xmean);
    xk = xk(2:end,:) + repmat(p.xmean',[N 1]); xk = xk'; 
elseif strcmp(select_model,'SINDYc')
    Ns = size(x,1);
    xk = zeros(Ns,N+1); xk(:,1) = x;
    for ct=1:N
        % Obtain plant state at next prediction step.
        xk(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xk(:,ct),u(ct),p.dt,1,[],p);
    end
    xk = xk(:,2:N+1);
elseif strcmp(select_model,'NARX')    
    Hu = [u',0];
    Hx = zeros(2,length(Hu)); Hx(:,1) = x;
    [Us,Ui,Si] = preparets(p.net,con2seq(Hu),{},con2seq(Hx));
    xk = p.net(Us,Ui,Si);
    xk = cell2mat(xk); 
end


%% Inequality constraints calculation
c = zeros(N,1);
% Apply N population size constraints across prediction horizon, from time
% k+1 to k+N
uk = u(1);
for ct=1:N
    % -z + zMin < 0 % lower bound
    c(ct) = -xk(2,ct)+zMin; %c(2*ct-1)
    % z - zMax < 0 % upper bound
    %c(2*ct) = xk1(1)-zMax;
    % update  input for next step
    if ct<N
        uk = u(ct+1);
    end
end
%% No equality constraints
ceq = [];

