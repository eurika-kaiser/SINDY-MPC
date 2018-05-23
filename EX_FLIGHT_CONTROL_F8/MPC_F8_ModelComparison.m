% Automatic flight control of F8 aircraft using model predictive control
% and various models

% MPC execution time for sine3 models
% 83.5243, 188.2009, 2.1348e+03
% 79.5731, 179.1829, 2.1387e+03

clear all, close all, clc

SystemModel = 'F8';

figpath = ['../FIGURES/',SystemModel,'/'];mkdir(figpath)
datapath = ['../DATA/',SystemModel,'/'];mkdir(datapath)
addpath('../utils');

%% Load Models
Nvar = 3;
InputSignalTypeModel = 'sine3';             % Use actuation to collect training data
InputSignalTypeModel_Validation = 'sine2';  % Use actuation to collect validation data to test for generalization

% DelayDMDc
ModelTypeDMDc = 'DMDc'; % Model: DelayDMDc , DMDc
load(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelTypeDMDc,'_',InputSignalTypeModel,'.mat'])) 
Models.DelayDMDc = Model;

% SINDYc
load(fullfile(datapath,['EX_',SystemModel,'_SI_SINDYc','_',InputSignalTypeModel,'.mat'])) 
Models.SINDYc = Model;

% NARX
NARX_SUBSTRACT_MEAN = 0;
load(fullfile(datapath,['EX_',SystemModel,'_SI_NARX','_',InputSignalTypeModel,'.mat'])) 
Models.NARX = Model;


%% Get training data
ENSEMBLE_DATA = 0;
ONLY_TRAINING_LENGTH = 1;
InputSignalType = InputSignalTypeModel;
Ndelay = Models.DelayDMDc.Ndelay;
getTrainingData

%% Get validation data: Model comparison
InputSignalType = InputSignalTypeModel_Validation;
Ndelay = Models.DelayDMDc.Ndelay;
tspanV =[tspan(end):dt:150];

% Forcing
if strcmp(InputSignalTypeModel_Validation,'sine2')
    forcing = @(x,t) [(.05* (sin(0.7*t).*sin(.1*t).*sin(.2*t).*sin(.05*t)) )];
elseif strcmp(InputSignalTypeModel_Validation,'sine3')
    forcing = @(x,t) (.5*sin(5*t).*sin(.5*t)+0.1).^3;
end

% Reference
[tA,xA] = ode45(@(t,x)F8Sys(t,x,forcing(x,t)),tspanV,x(end,1:3),options);   % true model
u_valid = forcing(0,tA)';

% DelayDMDc / DMDc
if Ndelay == 1
    x0      = [x(end,1:Nvar)];
    Hunew   = [u(end),u_valid(1:end-1)];
    [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspanV,[x0-[Models.DelayDMDc.xmean]]');
elseif Ndelay > 1
    x0      = [x(end-Ndelay+1,1:Nvar),x(end,1:Nvar)];
    Hunew   = [ u(end-Ndelay+1:end),u_valid(1:end-Ndelay);
        u(end),u_valid(1:end-1)];
    [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspanV,[x0-[Models.DelayDMDc.xmean]]');
    xB = xB(:,Nvar+1:2*Nvar); xB = xB + repmat(xmean,[size(xB,1) 1]);
end
xB = xB + repmat(Models.DelayDMDc.xmean,[length(tB) 1]);

% SINDYc
[tC,xC]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Models.SINDYc.Xi(:,1:Nvar),Models.SINDYc.polyorder,Models.SINDYc.usesine),tspanV,x(end,1:Nvar),options);  % approximate

% NARX
Udummy = [u(end),u_valid(1:end)];
Hdummy = zeros(Nvar,size(Udummy,2));
if NARX_SUBSTRACT_MEAN == 1
    NARX_xmean = Models.NARX.xmean;
else
    NARX_xmean = zeros(size(Models.NARX.xmean'));
end

Hdummy(:,1) = x(end,1:Nvar)'-NARX_xmean';
[Us,Ui,Si] = preparets(Models.NARX.net,con2seq(Udummy),{},con2seq(Hdummy));
xD = Models.NARX.net(Us,Ui,Si); % Predict on validation data
xD = cell2mat(xD)'; xD = [x0;xD(1:end-1,:)]; 
xD(2:end,:) = xD(2:end,:) + repmat(NARX_xmean,[size(xD,1)-1 1]);

%% Show results
clear ph
h = figure;
subplot(3,1,1), box on, hold on
plot([tB(1) tB(1)],[-100 100],'--k')
plot(t,x(:,1),'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,1),'k','LineWidth',1.5);
plot(tB,xB(:,1),'r-','LineWidth',1.5);
plot(tC,xC(:,1),'g--','LineWidth',1.5);
plot(tA(1:end),xD(:,1),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
grid on
ylim([-0.4 0.3])
ylabel('x_1','FontSize',13)
set(gca,'FontSize',13)

subplot(3,1,2), box on, hold on
plot([tB(1) tB(1)],[-100 100],'--k')
plot(t,x(:,2),'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,2),'k','LineWidth',1.5);
plot(tB,xB(:,2),'r-','LineWidth',1.5);
plot(tC,xC(:,2),'g--','LineWidth',1.5);
plot(tA(1:end),xD(:,2),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
ylim([-4 0.1])
set(gca,'FontSize',13)
ylabel('x_2','FontSize',13)

subplot(3,1,3), box on, hold on
plot([tB(1) tB(1)],[-100 100],'--k')
ph(1) = plot(t,x(:,3),'Color',[.4 .4 .4],'LineWidth',1.5);
ph(2) = plot(tA,xA(:,3),'k','LineWidth',1.5);
ph(3) = plot(tB,xB(:,3),'r-','LineWidth',1.5);
ph(4) = plot(tC,xC(:,3),'g--','LineWidth',1.5);
ph(5) = plot(tA(1:end),xD(:,3),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
l1=legend(ph,'Training','Validation',ModelTypeDMDc,'SINDYc','NARX');
set(l1,'Location','NorthWest')
grid on
ylim([-0.9 0.5])
ylabel('x_3','FontSize',13)
set(gca,'FontSize',13)
xlabel('Time','FontSize',13)


set(h,'Units','Inches');
set(gcf,'Position',[1 1 6. 5.5])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
print(h,'-painters','-depsc2',  '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'.eps'],'-r0');


%% Actuation signal
t_valid = tspanV;
clear ph
figure,box on, hold on,
ccolors = get(gca,'colororder');
plot([tA(1),tA(1)],[-15 260],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
plot(t,u,'-k','LineWidth',1);
plot(t_valid,u_valid,'-k','LineWidth',1);
grid off
ylim([min([u u_valid])+0.05*min([u u_valid]) max([u u_valid])+0.05*max([u u_valid])])
xlim([0 t_valid(end)])
xlabel('Time')
ylabel('Input')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_Actuation.eps']);

%% True SINDYc model parameters
Xi0 = zeros(size(Models.SINDYc.Xi));
Xi0([2,4,5,6,8,10,16,18,19,25,35],1) = [-0.877,1,-0.215,0.47,-0.088,-0.019,3.846,-1,0.28,0.47,0.63];
Xi0([4],2) = 1;
Xi0([2,4,5,6,16,19,25,35],3) = [-4.208,-0.396,-20.967,-0.47,-3.564,6.265,46,61.4];

%% MPC
% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
Duration = 6;                   % Run for 'Duration' time units
Ton = 0;                        % Time units when control is turned on
getMPCparams
          
x0n=[0.1 0 0]'; %xA(end,:)';    % Initial condition
Tcontrol = tA(end);             % Time offset to combine training, prediction and control phase

% Parameters Models
ModelCollection = {ModelTypeDMDc, 'SINDYc', 'NARX'}; % 51.8924, 46.0083, 1.4662e+03 execution times
Nmodels = length(ModelCollection);

clear Results
Results(1:Nmodels) = struct('x',[], 'u', [], 't', [], 'xref', [], 'J', [], 'elapsed_time', []);

Ts = Models.SINDYc.dt;

%%
for iM = 1:Nmodels
    select_model = ModelCollection{iM}; % DelayDMDc; DMDc; NARX ; SINDYc
    
    % Prepare variables
    Nt = (Duration/Ts)+1;
    uopt0    = 0;
    xhat     = x0n;
    uopt     = uopt0.*ones(Nu,1);
    xHistory = zeros(Nvar,Nt); xHistory(:,1) = xhat;
    uHistory = zeros(1,Nt); uHistory(1)   = uopt(1);
    tHistory = zeros(1,Nt); tHistory(1)   = Tcontrol;
    rHistory = zeros(1,Nt);
    
    % Parameters Model
    switch select_model
        case 'DelayDMDc'
            pest.dt     = Models.DelayDMDc.dt;
            pest.sys    = Models.DelayDMDc.sys;
            pest.xmean  = xmean';
            pest.udelay = zeros(1,N);
            pest.xdelay = [x(end-Ndelay+1,1:2)]';
            pest.Nxlim  = size(x,1);
        case 'DMDc'
            pest.dt = Models.DelayDMDc.dt;
            pest.sys = Models.DelayDMDc.sys;
            pest.xmean = Models.DelayDMDc.xmean';
            pest.Nxlim = size(x,1);  
        case 'SINDYc'
            pest.ahat = Models.SINDYc.Xi(:,1:Nvar); % Replace with Xi0 to test true model
            pest.polyorder = Models.SINDYc.polyorder;
            pest.usesine = Models.SINDYc.usesine;
            pest.dt = Models.SINDYc.dt;
        case 'NARX'
            pest.net = Models.NARX.net;
            pest.Ndelay = 1;
            pest.udelay = zeros(1,Ndelay);
            pest.xdelay = zeros(2,length(pest.udelay));
            pest.xdelay(:,1:Ndelay) = [x(end-Ndelay+1:end,1:2)]';
            pest.Nxlim = size(x,1);
    end
    
    
    % Start simulation
    fprintf('Simulation started.  It might take a while...\n')
    tic
    for ct = 1:(Duration/Ts)
        
        % Set references over prediction horizon
        tref = (ct:ct+N-1).*Ts;
        xref = [xrefFUN(tref); zeros(1,N); zeros(1,N)];
        
        % NMPC with full-state feedback
        COSTFUN = @(u) ObjectiveFCN_models(u,xhat,N,Nu,xref,uHistory(:,ct),pest,diag(Q),R,Ru,select_model);
        CONSFUN = @(u) ConstraintFCN_models(u,uHistory(:,ct),xhat,N,LBo,UBo,LBdu,UBdu,pest,select_model);
        uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
        
        % Run without constraints
%         uopt = fmincon(COSTFUN,uopt,[],[],[],[],[],[],[],options);        
        
        % Integrate system
        xhat = rk4u(@F8Sys,xhat,uopt(1),Ts/10,10,[],0); 
        xHistory(:,ct+1) = xhat;
        uHistory(:,ct+1) = uopt(1);
        tHistory(:,ct+1) = ct*Ts+Tcontrol;
        rHistory(:,ct+1) = xref(1,2);
       
    end
    tElapsed = toc    
    fprintf('Simulation finished!\n')
    
    % Collect results
    Results(iM).eval_time = tElapsed;
    Results(iM).xref = rHistory;
    Results(iM).x = xHistory;
    Results(iM).u = uHistory;
    Results(iM).t = tHistory;
    Results(iM).J = evalObjectiveFCN(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);
    Results(iM).elapsed_time = tElapsed;

    VIZ_SI_Validation_MPC
    
    if iM==Nmodels % Plot ensemble training for NN / NARX
       VIZ_SI_Validation_MPC_ensemble(ct,t_valid,u_valid,xHistory,tHistory,uHistory,Results,t,tA,xA,xB,xC,xD,select_model,SystemModel,figpath,N,InputSignalTypeModel,Nmodels,ModelCollection)
    end
    
end

% MPC Execution times using N=10: 120.9128, 96.8967, 2.1755e+03

%% Unforced
[tU,xU] = ode45(@(t,x)F8Sys(t,x,0,[]),tHistory,xA(end,1:3),options);   % true model

%% Show cost
clear ph
figure, hold on, box on
ph(1) = plot(Results(1).t,cumsum(Results(1).J),'-r','LineWidth',2);
ph(3) = plot(Results(1).t,cumsum(Results(3).J),'-.','Color',[0.7,0.7,1],'LineWidth',2);
ph(2) = plot(Results(1).t,cumsum(Results(2).J),'--g','LineWidth',2);

xlim([Results(1).t(1) Results(1).t(end)])
ylabel('Cost','FontSize',14)
xlabel('Time','FontSize',14)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'.eps']);
    

%% Show control stage separate
clear ph
figure, hold on, box on
ph(4) = plot(Results(1).t(2:end),(Results(1).xref(2:end)),'-k','LineWidth',2);
ph(1) = plot(Results(1).t,(Results(1).x(1,:)),'-r','LineWidth',2);
ph(3) = plot(Results(1).t,(Results(3).x(1,:)),'-.','Color',[0.7,0.7,1],'LineWidth',2);
ph(2) = plot(Results(1).t,(Results(2).x(1,:)),'--g','LineWidth',2);

l1=legend(ph,ModelCollection,'Ref');
set(l1,'Location','NorthEast')
xlim([Results(1).t(1) Results(1).t(end)])
ylim([-0.25 .25])
ylabel('AoA','FontSize',14)
xlabel('Time','FontSize',14)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_AoA','_N_',num2str(N),'.eps']);

clear ph
figure, hold on, box on
ph(1) = plot(Results(1).t,(Results(1).u(:)),'-r','LineWidth',2);
ph(3) = plot(Results(1).t,(Results(3).u(:)),'-.','Color',[0.7,0.7,1],'LineWidth',2);
ph(2) = plot(Results(1).t,(Results(2).u(:)),'--g','LineWidth',2);

xlim([Results(1).t(1) Results(1).t(end)])
ylabel('Input','FontSize',14)
xlabel('Time','FontSize',14)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Input','_N_',num2str(N),'.eps']);


%% Save results
save(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_N_',num2str(N),'.mat']),'Results')
save(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_All','_N_',num2str(N),'.mat']))