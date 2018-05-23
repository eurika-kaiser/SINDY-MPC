% LOTKA-VOLTERRA system
% System identification: DelayDMDc

clear all, close all, clc
figpath = '../FIGURES/';
datapath = '../DATA/';
addpath('../utils');

%% Generate Data %{'sine2', 'chirp','prbs', 'sphs'}
InputSignalType = 'sphs'; %prbs; chirp; noise; sine2; sphs; mixed
ONLY_TRAINING_LENGTH = 1;
DERIV_NOISE = 0;
getTrainingData

%% SINDYc
% Parameters
ModelName = 'SINDYc';
polyorder = 3;
usesine = 0;
lambda = 0.5; 
trainSINDYc

%% Prediction  over training phase
if any(strcmp(InputSignalType,{'sine2', 'chirp','prbs', 'sphs'})==1)
    [tSINDYc,xSINDYc]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xi(:,1:2),polyorder,usesine),tspan,x0,options);  % approximate
else
    p.ahat = Xi(:,1:2);
    p.polyorder = polyorder;
    p.usesine = usesine;
    p.dt = dt;
    [N,Ns] = size(x);
    xSINDYc = zeros(Ns,N); xSINDYc(:,1) = x0';
    for ct=1:N-1
        xSINDYc(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xSINDYc(:,ct),u(ct),dt,1,[],p);
    end
    xSINDYc = xSINDYc';
end

%% PhasePlots

clear ph
figure, hold on, box on
ccolors = get(gca,'colororder');
plot(x(:,1), x(:,2),'-k','LineWidth',1);
axis equal, axis tight, 
set(gca,'xtick', [], 'ytick', [])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_',InputSignalType,'_PhasePlot','_TRUTH.eps']);

clear ph
ccolors = [1 0 0; 0 1 0; 0.7,0.7,1];
figure, hold on, box on
plot(x(:,1), x(:,2),'-','Color','k','LineWidth',1.2);
plot(xSINDYc(:,1), xSINDYc(:,2),'--','Color',ccolors(2,:),'LineWidth',1.2);
axis equal, axis tight, %axis off
set(gca,'xtick', [], 'ytick', [])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 150])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_',InputSignalType,'_PhasePlot','_MODEL.eps']);

clear ph
figure, hold on, box on
ccolors = get(gca,'colororder');
iM = 1;
plot(tspan, u,'-k','LineWidth',1);
axis equal, axis tight, %axis off
set(gca,'xtick', [], 'ytick', [])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 50])
set(gcf,'PaperPositionMode','auto'),
print('-depsc2', '-painters', '-loose', '-cmyk', [figpath,'EX_LOTKA_',InputSignalType,'_Actuation','_TRUTH.eps']);


%% Show validation
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ph(1) = plot(tspan,x(:,1),'-','Color',ccolors(1,:),'LineWidth',1); hold on
ph(2) = plot(tspan,x(:,2),'-','Color',ccolors(2,:),'LineWidth',1);
ph(3) = plot(tspan,xSINDYc(:,1),'--','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
ph(4) = plot(tspan,xSINDYc(:,2),'--','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
xlim([0 100]), ylim([0 max(x(:))])
xlabel('Time')
ylabel('Population size')
legend(ph([1,3]),'True',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.eps']);

%% Prediction
% Reference
tspanV   = [100:dt:200];
xA      = xv;
tA      = tv;

% Model
if any(strcmp(InputSignalType,{'sine2', 'chirp','prbs', 'sphs'})==1)
    [tB,xB]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xi(:,1:2),polyorder,usesine),tspanV,x(end,:),options);  % approximate
    xB = xB(2:end,:);
    tB = tB(2:end); % rm IC
else
    [N,Ns] = size(xA);
    xB = zeros(Ns,N); xB(:,1) = x(end,:)';
    for ct=1:N
        xB(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xB(:,ct),uv(ct),dt,1,[],p);
    end
    xB = xB(:,1:N+1)';
    tB = tspanV(1:end);
end
%% Show training and prediction
VIZ_SI_Validation

%% Save Data
Model.name = 'SINDYc';
Model.polyorder = polyorder;
Model.usesine = usesine;
Model.Xi = Xi;
Model.dt = dt;
save(fullfile(datapath,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')