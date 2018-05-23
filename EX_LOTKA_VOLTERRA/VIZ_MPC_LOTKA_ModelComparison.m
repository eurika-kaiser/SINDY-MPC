% LOTKA-VOLTERRA system

clear all, close all, clc
figpath = '../FIGURES/NOT_USED/'; mkdir(figpath)
datapath = '../DATA/';
addpath('../utils');

%% Load Models
InputSignalTypeModel = 'sphs'; % prbs; chirp; noise; sine2; sphs; mixed

% DelayDMDc
ModelTypeDMDc = 'DMDc';
load(fullfile(datapath,['EX_LOTKA_SI_',ModelTypeDMDc,'_',InputSignalTypeModel,'.mat'])) % Model: DelayDMDc , DMDc
Models.DelayDMDc = Model;

% SINDYc
load(fullfile(datapath,['EX_LOTKA_SI_SINDYc','_',InputSignalTypeModel,'.mat'])) % Model
Models.SINDYc = Model;

% NARX
load(fullfile(datapath,['EX_LOTKA_SI_NARX','_',InputSignalTypeModel,'.mat'])) % Model
Models.NARX = Model;

%% Get training data
InputSignalType = InputSignalTypeModel;% prbs; chirp; noise; sine2; sphs; mixed
Ndelay = Models.DelayDMDc.Ndelay;
ONLY_TRAINING_LENGTH = 1;
getTrainingData

%% Get validation data
InputSignalType = 'sine2';% prbs; chirp; noise; sine2; sphs; mixed
Ndelay = Models.DelayDMDc.Ndelay;
getValidationData


%% FIGURE 1:  // Model comparison + Validation

% Reference
[tA,xA] = ode45(@(t,x)lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspanV,x(end,1:2),options);   % true model

% DelayDMDc / DMDc
% Model
if Ndelay == 1
    x0      = [x(end,1:2)];
    Hunew   = [u(end),u_valid(1:end-1)];
    [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspanV,[x0-[xmean]]');
elseif Ndelay > 1
    x0      = [x(end-Ndelay+1,1:2),x(end,1:2)];
    Hunew   = [ u(end-Ndelay+1:end),u_valid(1:end-Ndelay);
        u(end),u_valid(1:end-1)];
    [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspanV,[x0-[xmean]]');
    xB = xB(:,3:4); xB = xB + repmat(xmean,[size(xB,1) 1]);
end
xB = xB + repmat(xmean,[length(tB) 1]);

% SINDYc
[tC,xC]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Models.SINDYc.Xi(:,1:2),Models.SINDYc.polyorder,Models.SINDYc.usesine),tspanV,x(end,1:2),options);  % approximate

% NARX
Udummy = [u(end),u_valid(1:end)];
Hdummy = zeros(2,size(Udummy,2));
Hdummy(:,1) = x(end,1:2)';
[Us,Ui,Si] = preparets(Models.NARX.net,con2seq(Udummy),{},con2seq(Hdummy));
xD = Models.NARX.net(Us,Ui,Si); % Predict on validation data
xD = cell2mat(xD)'; xD = [x0;xD(1:end-1,:)];

%%
clear ph
h = figure;
subplot(2,1,1), box on, hold on
plot([tB(1) tB(1)],[0 150],'--k')
plot(t,x(:,1),'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,1),'k','LineWidth',1.5);
plot(tB,xB(:,1),'r-','LineWidth',1.5);
plot(tC,xC(:,1),'g--','LineWidth',1.5);
plot(tA(1:end),xD(:,1),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
grid on
ylim([0 150])
% xlabel('Time','FontSize',13)
ylabel('Prey, x_1','FontSize',13)
set(gca,'FontSize',13)
subplot(2,1,2), box on, hold on
plot([tB(1) tB(1)],[0 60],'--k')
ph(1) = plot(t,x(:,2),'Color',[.4 .4 .4],'LineWidth',1.5);
ph(2) = plot(tA,xA(:,2),'k','LineWidth',1.5);
ph(3) = plot(tB,xB(:,2),'r-','LineWidth',1.5);
ph(4) = plot(tC,xC(:,2),'g--','LineWidth',1.5);
% ph(5) = plot(tB,xD(:,2),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
ph(5) = plot(tA(1:end),xD(:,2),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
l1=legend(ph,'Training','Validation',ModelTypeDMDc,'SINDYc','NARX');
set(l1,'Location','NorthWest')
grid on
ylim([0 60])
ylabel('Predator, x_2','FontSize',13)
set(gca,'FontSize',13)
xlabel('Time','FontSize',13)
set(gca,'FontSize',13)

set(h,'Units','Inches');
set(gcf,'Position',[1 1 6. 5.5])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
print(h,'-painters','-depsc2',  '-loose','-cmyk', [figpath,'EX_LOTKA_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'.eps'],'-r0');


%% Actuation signal
clear ph
figure,box on, hold on,
ccolors = get(gca,'colororder');
plot([tA(1),tA(1)],[-15 260],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
plot(t,u,'-k','LineWidth',1);
plot(t_valid,u_valid,'-k','LineWidth',1);
grid off
ylim([min([u u_valid])+0.05*min([u u_valid]) max([u u_valid])+0.05*max([u u_valid])])
xlim([0 200])
xlabel('Time')
ylabel('Input')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_LOTKA_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_Actuation.eps']);

%% MPC
% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
N  = 10;                        % Control / prediction horizon (number of iterations)
Duration = 100;                 % Run for 'Duration' time units
Ton = 0;                        % Time units when control is turned on
Nvar = 2;
getMPCparams
                
x0n=xA(end,:)';                 % Initial condition
xref1 = [g/d;a/b];              % critical point
Tcontrol = tA(end);             % Time offset to combine training, prediction and control phase

% Parameters Models
ModelCollection = {ModelTypeDMDc, 'SINDYc', 'NARX'};
Nmodels = length(ModelCollection);

% Parameters True Model
p.a = a;
p.b = b;
p.d = d;
p.g = g;

%% Load Results
clear Results
load(fullfile(datapath,['MPC_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'.mat']),'Results')


%% Unforced
[tU,xU] = ode45(@(t,x)lotkacontrol(t,x,0,a,b,d,g),Results(1).t,xA(end,1:2),options);   % true model

%% Show results
MODEL_COLOR = 1;
Ts = Models.DelayDMDc.dt;
ct = (Duration/Ts);
for iM = 1:Nmodels
    select_model = ModelCollection{iM};
    xHistory = Results(iM).x;
    uHistory = Results(iM).u;
    tHistory = Results(iM).t;
    VIZ_SI_Validation_MPC
end
     

%% Get cost // evaluate performance
% Data2 = Datal
Costs = zeros(length(Results(1).u),3,Nmodels);

for iM = 1:Nmodels
    [J] = evalObjectiveFCN(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);
    Costs(:,1,iM) = J; %Costs(:,2,iM) = Js; Costs(:,3,iM) = Ju;
end

%% Show cost
clear ph
figure, hold on, box on
ph(1) = plot(Results(1).t,cumsum(Results(1).J),'-r','LineWidth',2);
ph(2) = plot(Results(1).t,cumsum(Results(2).J),'--g','LineWidth',2);
ph(3) = plot(Results(1).t,cumsum(Results(3).J),'-.','Color',[0.7,0.7,1],'LineWidth',2);

xlim([200 xlimval(2)])
l1=legend(ph,ModelCollection);
set(l1,'Location','SouthEast')
% ylim([0 5.5*10^4])
ylabel('Cost','FontSize',14)
xlabel('Time','FontSize',14)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_LOTKA_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_Cost.eps']);

clear ph
figure, hold on, box on
ph(1) = plot(Results(1).t,cumsum(Costs(:,1,1)),'-r','LineWidth',2);
ph(2) = plot(Results(1).t,cumsum(Costs(:,1,2)),'--g','LineWidth',2);
ph(3) = plot(Results(1).t,cumsum(Costs(:,1,3)),'-.','Color',[0.7,0.7,1],'LineWidth',2);

xlim([200 xlimval(2)])
l1=legend(ph,ModelCollection);
set(l1,'Location','SouthEast')
% ylim([0 5.5*10^4])
ylabel('Cost','FontSize',14)
xlabel('Time','FontSize',14)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
%%
iJ = 3;
clear ph
figure, hold on, box on
ph(1) = plot(Results(1).t,cumsum(Costs(:,iJ,1)),'-r','LineWidth',2);
ph(2) = plot(Results(1).t,cumsum(Costs(:,iJ,2)),'--g','LineWidth',2);
ph(3) = plot(Results(1).t,cumsum(Costs(:,iJ,3)),'-.','Color',[0.7,0.7,1],'LineWidth',2);

l1=legend(ph,ModelCollection);
set(l1,'Location','NorthEast')
ylabel('Cost','FontSize',14)
xlabel('Time','FontSize',14)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
    
