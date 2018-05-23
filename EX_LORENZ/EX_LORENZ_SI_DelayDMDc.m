% LORENZ system
% System identification: DelayDMDc

clear all, close all, clc
figpath = '../FIGURES/LORENZ/';
datapath = '../DATA/LORENZ/';
addpath('../utils');

SystemModel = 'LORENZ';

%% Generate Data
InputSignalType = 'sphs';%prbs; chirp; noise; sine2; sphs; mixed
Ndelay = 1;
ONLY_TRAINING_LENGTH = 1;
getTrainingData

Nt = length(tspan)-1;

%% DMDc: B = unknown  and with time delay coordinates
ModelNumber = 2; % or 2
xrefs = [xref1,xref2];

if Ndelay == 1
    ModelName = 'DMDc';
elseif Ndelay>1
    ModelName = 'DelayDMDc';
end

Hu = getHankelMatrix_MV(u',1);
for i = 1:ModelNumber
    xmean{i} = xrefs(:,i)'; %mean(x);
    X   = x - repmat(xmean{i},[T 1]);
    Hx  = getHankelMatrix_MV(X,1);
    numOutputs = size(Hx,1); numInputs = size(Hu,1); numVar = 3;
    r1 = size(Hx,1); r2 = size(Hx,1);
    [sysmodel_DMDc{i},U,Up] = DelayDMDc_MV(Hx,Hu,size(Hx,1),size(Hx,1),dt,size(Hx,1),size(Hu,1),2);
end
%% Prediction over training phase
for i = 1:ModelNumber
    [xDMDc{i},~] = lsim(sysmodel_DMDc{i},Hu',tspan(1:end),x(1,:)'-xmean{i}');
    xDMDc{i} = xDMDc{i} + repmat(xmean{i},[length(tspan) 1]);
end

%% Show validation
for i = 1:ModelNumber
    clear ph
    figure,box on,
    ccolors = get(gca,'colororder');
    ph(1) = plot(tspan,x(:,1),'-','Color',ccolors(1,:),'LineWidth',1); hold on
    ph(2) = plot(tspan,x(:,2),'-','Color',ccolors(2,:),'LineWidth',1);
    ph(3) = plot(tspan(Ndelay:Nt+Ndelay),xDMDc{i}(:,1),'--','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
    ph(4) = plot(tspan(Ndelay:Nt+Ndelay),xDMDc{i}(:,2),'--','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
    xlim([0 (length(tspan)-1)*dt]), ylim([-25 50])
    xlabel('Time')
    ylabel('Population size')
    % legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
    legend(ph([1,3]),'True',ModelName)
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_M',i,'.eps']);
end
%% Prediction
% Reference
tspanV   = [10:dt:20];
xA      = xv;
tA      = tv;

% Model
for i = 1:ModelNumber
    if Ndelay == 1
        x0      = [x(end,1:3)];
        Hunew   = [u(end),uv(1:end)];
        [xBm{i},tB] = lsim(sysmodel_DMDc{i},Hunew,tspanV,[x0-[xmean{i}]]');
    elseif Ndelay > 1
        x0      = [x(end-Ndelay+1,1:3),x(end,1:3)];
        Hunew   = [ u(end-Ndelay+1:end),uv(1:end-Ndelay);
            u(end),uv(1:end-1)];
        [xBm,tB] = lsim(sysmodel_DMDc,Hunew,tspanV,[x0-[xmean]]');
        xBm = xBm(:,4:6); xBm = xBm + repmat(xmean,[size(xBm,1) 1]);
    end
    
    xBm{i} = xBm{i} + repmat(xmean{i},[length(tB) 1]);
end
%% Show training and prediction
xB = xBm{1};
VIZ_SI_Validation

xB = xBm{2};
VIZ_SI_Validation

%% Save Data
Model.name = 'DelayDMDc';
Model.sys = sysmodel_DMDc;
Model.Ndelay = Ndelay;
Model.xmean = xmean;
Model.xrefs = xrefs;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')