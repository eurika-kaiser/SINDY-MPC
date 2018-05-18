function [c, ceq] = ConstraintFCN_models(u,uold,x,N,LBo,UBo,LBdu,UBdu,p,select_model)
%% Constraint function of nonlinear MPC for F8 system
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon length
%   uold:   latest applied control input
%   LBo:    Lower bound of output x
%   UBo:    Upper bound of output x
%   LBdu:   Lower bound for input difference uk - uk-1
%   UBdu:   Upper bound for input difference uk - uk-1
%   p:      Parameters for model
%   select_model: Selects model future-state prediction
%
% Output:
%   c:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints  (empty)
%

Nvar = length(x);
%% Nonlinear MPC design parameters
% Ensure that all cell populations are positive
zMin = LBo; 

%% Integrate system
if strcmp(select_model,'DelayDMDc')
    Uinput = getHankelMatrix_MV([p.udelay(2:end)';u],p.Ndelay)';
    [xk,~] = lsim(p.sys,Uinput,[0:N-1].*p.dt,p.xdelay);
    xk = xk(:,end-Nvar+1:end);
    xk = xk + repmat(p.xmean',[N 1]); xk = xk';    
elseif strcmp(select_model,'DMDc')
    [xk,~] = lsim(p.sys,[u',0],[0:N].*p.dt,x-p.xmean);
    xk = xk(2:end,:) + repmat(p.xmean',[N 1]); xk = xk'; 
elseif strcmp(select_model,'eDMDc')
    Y0 = poolData((x-p.xmean)',Nvar,p.polyorder,p.usesine)';
    [xk,~] = lsim(p.sys,[u' 0],[0:N].*p.dt,Y0(2:end));
    xk = xk(2:end,1:Nvar) + repmat(p.xmean',[N 1]); xk = xk';       
elseif strcmp(select_model,'SINDYc')
    Ns = size(x,1);
    xk = zeros(Ns,N+1); xk(:,1) = x;
    for ct=1:N
        % Obtain plant state at next prediction step.
        xk(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xk(:,ct),u(ct),p.dt,1,[],p);
    end
    xk = xk(:,2:N+1);
elseif strcmp(select_model,'partialSINDYc')
    Ns = size(x,1);
    xk = zeros(Ns,N+1); xk(:,1) = x;
    for ct=1:N
        % Obtain plant state at next prediction step.
        xk_tmp = rk4u(@sparseGalerkinControl_Discrete,xk(p.SelectVars,ct),u(ct),p.dt,1,[],p);
        xk(p.SelectVars,ct+1) = xk_tmp;
    end
    xk = xk(:,2:N+1);     
elseif strcmp(select_model,'NARX')    
    Hu = [u',0];
    Hx = zeros(Nvar,length(Hu)); Hx(:,1) = x;
    if p.TRANSFORM_LOG == 1
        Hx = log(Hx);
    end
    [Us,Ui,Si] = preparets(p.net,con2seq(Hu),{},con2seq(Hx));
    xk = p.net(Us,Ui,Si);
    xk = cell2mat(xk); 
    if p.TRANSFORM_LOG == 1
        xk = exp(xk);
    end
end

    

%% Inequality constraints calculation
c = zeros(N,1);
% Apply N population size constraints across prediction horizon, from time
% k+1 to k+N

for ct=1:N
    % -z + zMin < 0 % lower bound
    c(ct) = -xk(1,ct)+zMin;

end

%% No equality constraints
ceq = [];

