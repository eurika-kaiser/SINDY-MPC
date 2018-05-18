function J = ObjectiveFCN_models(u,x,N,xref,u0,p,Q,R,Ru, select_model)
%% Cost function of nonlinear MPC for a dynamical system
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

Nvar = length(x);
%% Nonlinear MPC design parameters
% Q = diag([1,1,1]);
% R = 0.01;

%% Integrate system
if strcmp(select_model,'DelayDMDc')
%     Hunew   = [ u(end-Ndelay+1:end),unew(1:end-Ndelay);
%     u(end),unew(1:end-1)];

%     [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspan,[x0-[xmean,xmean]]');
%     xB = xB(:,3:4); xB = xB + repmat(xmean,[size(xB,1) 1]);

%     Ns = size(x,1);
%     xk = zeros(2*Ns,N+1); xk(:,1) = [p.xdelay; x]-[p.xmean; p.xmean];
%     for ct=1:N
%         % Obtain plant state at next prediction step.
%         xk(:,ct+1) = p.sys.A*xk(:,ct) + p.sys.B*[p.udelay(ct);u(ct)];
%     end
%     xk = xk(3:4,:);
%     xk = xk + repmat(p.xmean,[1 N+1]);
%     xk = xk(:,2:N+1);
    
    [xk,~] = lsim(p.sys,[p.udelay(1:N);u'],[0:N-1].*p.dt,[p.xdelay(:,1); x]-[p.xmean; p.xmean]);
    xk = xk(:,Nvar+1:2*Nvar);
    xk = xk + repmat(p.xmean',[N 1]); xk = xk';
elseif strcmp(select_model,'DMDc')
    [xk,~] = lsim(p.sys,[u' 0],[0:N].*p.dt,x-p.xmean);
    xk = xk(2:end,:) + repmat(p.xmean',[N 1]); xk = xk';
    
%     Ns = size(x,1); % identical
%     xk1 = zeros(Ns,N+1); xk1(:,1) = x-p.xmean;
%     for ct=1:N
%         % Obtain plant state at next prediction step.
%         xk1(:,ct+1) = p.sys.A*xk1(:,ct) + p.sys.B*u(ct);
%     end
%     xk1 = xk1(:,2:N+1)  + repmat(p.xmean,[1 N]);
elseif strcmp(select_model,'SINDYc')
    Ns = size(x,1);
    xk = zeros(Ns,N+1); xk(:,1) = x;
    for ct=1:N
        % Obtain plant state at next prediction step.
        xk(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xk(:,ct),u(ct),p.dt,1,[],p);
    end
    xk = xk(:,2:N+1);
elseif strcmp(select_model,'NARX')
%     Udummy = [p.udelay,u'];
%     Hdummy = zeros(2,size(Udummy,2));
%     Hdummy(:,1:p.Ndelay) = p.xdelay- repmat(p.xmean,[1 p.Ndelay]); %length(p.udelay)
%     [Us,Ui,Si] = preparets(p.net,con2seq(Udummy),{},con2seq(Hdummy));
%     xk = p.net(Us,Ui,Si);
%     xk = cell2mat(xk); xk = xk + repmat(p.xmean,[1 size(xk,2)]);

%     Hu = [u'-p.umean,0];
    Hu = [u',0];
    Hx = zeros(Nvar,length(Hu)); Hx(:,1) = x;%-p.xmean;
    [Us,Ui,Si] = preparets(p.net,con2seq(Hu),{},con2seq(Hx));
    xk = p.net(Us,Ui,Si);
    xk = cell2mat(xk); xk = xk;% + repmat(p.xmean,[1 size(xk,2)]);
%     xk = [x,xk];
end

% figure(1);
% plot(u)
% pause(0.01)
%% Cost Calculation
% Set initial plant states, controller output and cost.
uk = u(1);
J = 0;
% Loop through each prediction step.
for ct=1:N
    % Obtain plant state at next prediction step.
    xk1 = xk(:,ct);
    
    % accumulate state tracking cost from x(k+1) to x(k+N).
    J = J + (xk1-xref)'*Q*(xk1-xref);
    % accumulate MV rate of change cost from u(k) to u(k+N-1).
    if ct==1
        J = J + (uk-u0)'*R*(uk-u0) + uk'*Ru*uk;
    else
        J = J + (uk-u(ct-1))'*R*(uk-u(ct-1)) + uk'*Ru*uk;
    end
    % Update uk for the next prediction step.
    if ct<N
        uk = u(ct+1);
    end
end
% J
