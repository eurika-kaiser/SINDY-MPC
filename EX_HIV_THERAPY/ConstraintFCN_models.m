function [c, ceq] = ConstraintFCN_models(u,uold,x,N,LBo,UBo,LBdu,UBdu,p,select_model)
%% Constraint function of nonlinear MPC for F8 system
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
Nvar = length(x);
%% Nonlinear MPC design parameters
% range of angle of attack
zMin = LBo;
zMax = UBo;

%% Integrate system
if strcmp(select_model,'DelayDMDc')
    Uinput = getHankelMatrix_MV([p.udelay(2:end)';u],p.Ndelay)';
    [xk,~] = lsim(p.sys,Uinput,[0:N-1].*p.dt,p.xdelay);
    xk = xk(:,end-Nvar+1:end);
    xk = xk + repmat(p.xmean',[N 1]); xk = xk';    
%      [xk,~] = lsim(p.sys,[p.udelay(1:N);u'],[0:N-1].*p.dt,[p.xdelay(:,1); x]-[p.xmean; p.xmean]);
%     xk = xk(:,Nvar+1:2*Nvar);
%     xk = xk + repmat(p.xmean',[N 1]); xk = xk';
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
    Hu = [u',0];%-p.umean,0];
    Hx = zeros(Nvar,length(Hu)); Hx(:,1) = x;%-p.xmean;
    if p.TRANSFORM_LOG == 1
        Hx = log(Hx);
    end
    [Us,Ui,Si] = preparets(p.net,con2seq(Hu),{},con2seq(Hx));
    xk = p.net(Us,Ui,Si);
    xk = cell2mat(xk); 
    if p.TRANSFORM_LOG == 1
        xk = exp(xk);
    end
    xk = xk;% + repmat(p.xmean,[1 size(xk,2)]);
end

    

%% Inequality constraints calculation
c = zeros(N,1);
% c = zeros(2*N,1);
% Apply N population size constraints across prediction horizon, from time
% k+1 to k+N
uk = u(1);
duk = u(1)-uold;
for ct=1:N
    % -z + zMin < 0 % lower bound
    c(ct) = -xk(1,ct)+zMin;
%     c(2*ct-1) = -xk(1,ct)+zMin; %c(ct)
%     c(2*ct-1) = -duk+LBdu; 
    % z - zMax < 0 % upper bound
%     c(2*ct) = xk(1,ct)-zMax;
%     c(2*ct) = duk-UBdu;
    % update  input for next step
    if ct<N
        uk = u(ct+1);
        duk = u(ct+1)-u(ct);
    end
end

% % cy = zeros(2*N,1); % upper AND lower bound
% cy = zeros(N,1); % upper OR lower bound
% for ct=1:N
%     % -z + zMin < 0 % lower bound
%     cy(ct) = -xk(1,ct)+zMin; %c(ct)
% %     cy(2*ct-1) = -xk(1,ct)+zMin; %c(ct)
%     % z - zMax < 0 % upper bound
% %     cy(2*ct) = xk(1,ct)-zMax;
% end
% 
% c = [c;cy];
%% No equality constraints
ceq = [];

