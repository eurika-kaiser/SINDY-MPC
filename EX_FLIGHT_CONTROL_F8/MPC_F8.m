% LORENZ system
%118.2675
%98.2860
%2.1048e+03

clear all, close all, clc

SystemModel = 'F8';

figpath = ['../../FIGURES/',SystemModel,'/'];mkdir(figpath)
datapath = ['../../DATA/',SystemModel,'/'];mkdir(datapath)
% datapath = ['/Users/ekaiser/Documents/Academia/Papers/KaKuBr_SINDYc-MPC/DATA/'];
addpath('../utils');

%% MPC


dt = 0.01;
Ts = 0.01;

% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
Duration = 6;                 % Run for 'Duration' time units
Ton = 0;                       % Time units when control is turned on
getMPCparams

x0n=[0.1 0 0]';                % Initial condition
Tcontrol = 0;                  % Time offset to combine training, prediction and control phase


clear Results
Results = struct('r', [], 'x',[], 'u', [], 't', [], 'xref', [], 'J', []);

%%
Nvar = 3;
% Prepare variables
Nt = (Duration/Ts)+1;
uopt0    = 0;
xhat     = x0n;
uopt     = uopt0.*ones(Nu,1);
xHistory = zeros(Nvar,Nt); xHistory(:,1) = xhat;
uHistory = zeros(1,Nt); uHistory(1)   = uopt(1);
tHistory = zeros(1,Nt); tHistory(1)   = Tcontrol;
rHistory = zeros(1,Nt);

% Start simulation
fprintf('Simulation started.  It might take a while...\n')
tic
for ct = 1:Duration/Ts %(Duration/Ts)
    % Set references
    tref = (ct:ct+N-1).*Ts;
    xref = [xrefFUN(tref); zeros(1,N); zeros(1,N)];
    %         figure(1);
    %         plot(tref,xref(1,:),'-k'),drawnow
    %         pause(0.1)
    
    % NMPC with full-state feedback
    COSTFUN = @(u) ObjectiveFCN(u,xhat,Ts,N,Nu,xref,uHistory(:,ct),[],diag(Q),R,Ru);
    CONSFUN = @(u) ConstraintFCN(u,uHistory(:,ct),xhat,Ts,N,LBo,UBo,LBdu,UBdu,[]);    
%     COSTFUN = @(u) ObjectiveFCN(u,xhat,Ts,N,Nu,xref,uHistory(:,ct+1),[],diag(Q),R,Ru);
%     CONSFUN = @(u) ConstraintFCN(u,uHistory(:,ct+1),xhat,Ts,N,LBo,UBo,LBdu,UBdu,[]);
    %         uopt = fmincon(COSTFUN,uopt,[],[],[],[],[],[],[],options);
    uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
    
    % Integrate system
    xhat = rk4u(@F8Sys,xhat,uopt(1),Ts/10,10,[],0); %10, 2
    xHistory(:,ct+1) = xhat;
    uHistory(:,ct+1) = uopt(1);
    tHistory(:,ct+1) = ct*Ts+Tcontrol;
    rHistory(:,ct+1) = xref(1,2);
    
end

fprintf('Simulation finished!\n')
%
tElapsed = toc

Results.eval_time = tElapsed;
Results.xref = rHistory;
Results.x = xHistory;
Results.u = uHistory;
Results.t = tHistory;
Results.J = evalObjectiveFCN(Results.u,Results.x,Results.xref,diag(Q),R,Ru);

figure,plot(xHistory'), hold on, plot(r,'-k'), legend('aoa','pa','pr','r')
figure,plot(uHistory)
figure,plot(Results.J)
% VIZ_SI_Validation_MPC

% [J,Js,Ju] = evalObjectiveFCN(Results.u,Results.x,Results.xref,diag(Q),R,Ru);

%%
% figure,hold on
% plot([1 length(Results(1).J); 1 length(Results(1).J); 1 length(Results(1).J)]',[xref1,xref1]','-k')
% plot(Results(1).x','-'), plot(Results(2).x','--')
%
% figure,plot(cumsum(Results(1).J'),'-'), hold on, plot(cumsum(Results(2).J'),'--')
% figure,plot(Results(1).u','-'), hold on, plot(Results(2).u','--')
