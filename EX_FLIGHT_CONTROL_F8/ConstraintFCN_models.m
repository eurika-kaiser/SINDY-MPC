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
% range of angle of attack
zMin = LBo(1);
zMax = UBo(1);

%% Integrate system
if strcmp(select_model,'DelayDMDc')
     [xk,~] = lsim(p.sys,[p.udelay(1:N);u'],[0:N-1].*p.dt,[p.xdelay(:,1); x]-[p.xmean; p.xmean]);
    xk = xk(:,Nvar+1:2*Nvar);
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
    Hx = zeros(Nvar,length(Hu)); Hx(:,1) = x;
    [Us,Ui,Si] = preparets(p.net,con2seq(Hu),{},con2seq(Hx));
    xk = p.net(Us,Ui,Si);
    xk = cell2mat(xk); 
end


%% Inequality constraints calculation
c = zeros(2*N,1);

% Apply N population size constraints across prediction horizon, from time
% k+1 to k+N
duk = u(1)-uold;
for ct=1:N
    c(2*ct-1) = -duk+LBdu; 
    c(2*ct) = duk-UBdu;
    if ct<N
        duk = u(ct+1)-u(ct);
    end
end

c1 = zeros(2*N,1);
c2 = zeros(2*N,1);
c3 = zeros(2*N,1);
for ct=1:N
    c1(2*ct-1) = -xk(1,ct)+zMin; 
    c2(2*ct-1) = -xk(2,ct)+LBo(2); 
    c3(2*ct-1) = -xk(3,ct)+LBo(3);
    c1(2*ct) = xk(1,ct)-zMax;
    c2(2*ct) = xk(2,ct)-UBo(2); 
    c3(2*ct) = xk(3,ct)-UBo(3);
end

c = [c;c1;c2;c3];

%% No equality constraints
ceq = [];

