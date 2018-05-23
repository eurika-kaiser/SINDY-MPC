% LOTKA-VOLTERRA system
% System identification: DelayDMDc

clear all, close all, clc
figpath = '../FIGURES/';
datapath = '../DATA/';
addpath('../utils');

%% Generate Data
InputSignalType = 'sine2';%prbs; chirp; noise; sine2; sphs; mixed
Ndelay = 1;
getTrainingData


%% DMDc: B = unknown  and with time delay coordinates
if Ndelay == 1
    ModelName = 'DMDc';
elseif Ndelay>1
    ModelName = 'DelayDMDc';
end

numOutputs = size(Hx,1); numInputs = size(Hu,1); numVar = 2;
r1 = size(Hx,1); r2 = size(Hx,1);
[sysmodel_DMDc,U,Up] = DelayDMDc_MV(Hx,Hu,size(Hx,1),size(Hx,1),dt,size(Hx,1),size(Hu,1),2);

%% Prediction over training phase
[xDMDc,~] = lsim(sysmodel_DMDc,Hu,tspan(1:Nt),Hx(:,1));
xDMDc = xDMDc(:,end-1:end);
xDMDc = xDMDc + repmat(xmean,[Nt 1]);


%% Show validation
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ph(1) = plot(tspan,x(:,1),'-','Color',ccolors(1,:),'LineWidth',1); hold on
ph(2) = plot(tspan,x(:,2),'-','Color',ccolors(2,:),'LineWidth',1);
ph(3) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,1),'--','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
ph(4) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,2),'--','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
xlim([0 100])
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
if Ndelay == 1
    x0      = [x(end,1:2)];
    Hunew   = [u(end),uv(1:end)];
    [xB,tB] = lsim(sysmodel_DMDc,Hunew,tspanV,[x0-[xmean]]');
elseif Ndelay > 1
    x0      = [x(end-Ndelay+1,1:2),x(end,1:2)];
    Hunew   = [ u(end-Ndelay+1:end),uv(1:end-Ndelay);
        u(end),uv(1:end-1)];
    [xB,tB] = lsim(sysmodel_DMDc,Hunew,tspanV,[x0-[xmean]]');
    xB = xB(:,3:4); xB = xB + repmat(xmean,[size(xB,1) 1]);
end

xB = xB + repmat(xmean,[length(tB) 1]);

%% Show training and prediction
VIZ_SI_Validation

%% Save Data
Model.name = 'DelayDMDc';
Model.sys = sysmodel_DMDc;
Model.Ndelay = Ndelay;
Model.dt = dt;
save(fullfile(datapath,['EX_LOTKA_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')