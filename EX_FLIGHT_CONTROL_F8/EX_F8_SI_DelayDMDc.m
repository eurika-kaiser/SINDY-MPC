% LOTKA-VOLTERRA system
% System identification: DelayDMDc

clear all, close all, clc
figpath = '../../FIGURES/F8/'; mkdir(figpath)
datapath = '../../DATA/F8/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'F8';
Nvar = 3;
%% Generate Data
InputSignalType = 'sine3';%prbs; chirp; noise; sine2; sine3; sphs; mixed
Ndelay = 1;
ENSEMBLE_DATA = 0;
ONLY_TRAINING_LENGTH = 1;
getTrainingData

Nt = length(tspan)-1;
xref = zeros(Nvar,1);
T = length(tspan);
%% DMDc: B = unknown  and with time delay coordinate

if Ndelay == 1
    ModelName = 'DMDc';
elseif Ndelay>1
    ModelName = 'DelayDMDc';
end

Hu = getHankelMatrix_MV(u',1);
xmean = xref'; %mean(x);
X   = x - repmat(xmean,[T 1]);
Hx  = getHankelMatrix_MV(X,1);
numOutputs = size(Hx,1); numInputs = size(Hu,1); numVar = 3;
r1 = size(Hx,1); r2 = size(Hx,1);
[sysmodel_DMDc,U,Up] = DelayDMDc_MV(Hx,Hu,size(Hx,1),size(Hx,1),dt,size(Hx,1),size(Hu,1),2);

%% Prediction over training phase
[xDMDc,~] = lsim(sysmodel_DMDc,Hu',tspan(1:end-1),x(1,:)'-xmean');
xDMDc = xDMDc + repmat(xmean,[length(tspan)-1 1]);

%% Show validation
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];ccolors(2,:)-[0.1 0.2 0.09];ccolors(3,:)-[0.1 0.2 0.09]];
for i = 1:Nvar
    ph(i) = plot(tspan,x(:,i),'-','Color',ccolors(i,:),'LineWidth',1); hold on
end
for i = 1:Nvar
    ph(Nvar+i) = plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
end

xlim([0 (length(tspan)-1)*dt]), ylim([-0.8 0.8])
xlabel('Time')
ylabel('Population size')
% legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
legend(ph([1,3]),'True',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);

%% Prediction
% Reference
tspanV   = [100:dt:200];
xA      = xv;
tA      = tv;

% Model

if Ndelay == 1
    x0      = [x(end,1:3)];
    Hunew   = [u(end),uv(1:end)];
    [xBm,tB] = lsim(sysmodel_DMDc,Hunew,tspanV,[x0-[xmean]]');
elseif Ndelay > 1
    x0      = [x(end-Ndelay+1,1:3),x(end,1:3)];
    Hunew   = [ u(end-Ndelay+1:end),uv(1:end-Ndelay);
        u(end),uv(1:end-1)];
    [xBm,tB] = lsim(sysmodel_DMDc,Hunew,tspanV,[x0-[xmean]]');
    xBm = xBm(:,4:6); xBm = xBm + repmat(xmean,[size(xBm,1) 1]);
end

xBm = xBm + repmat(xmean,[length(tB) 1]);

%% Show training and prediction
xB = xBm;
VIZ_SI_Validation

%% Save Data
Model.name = 'DelayDMDc';
Model.sys = sysmodel_DMDc;
Model.Ndelay = Ndelay;
Model.xmean = xmean;
Model.xrefs = xref;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')