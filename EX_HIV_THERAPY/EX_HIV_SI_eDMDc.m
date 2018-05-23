% HIV system
% System identification: extended DMD with control

clear all, close all, clc
figpath = '../FIGURES/HIV/'; mkdir(figpath)
datapath = '../DATA/HIV/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'HIV';

%% Generate Data
InputSignalType = 'prbs';
ONLY_TRAINING_LENGTH = 1;
DATA_ENSEMBLE = 0;

Nvar = 5;

ModelName = 'eDMDc';
getTrainingData


%% eDMDc: B = unknown  and with time delay coordinates
polyorder = 3;
usesine = 0;
Ndelay = 1;

% Construct data matrices
Hu = getHankelMatrix_MV(u,Ndelay);
xmean = xref;% zeros(1,Nvar); %mean(x); % In paper: zeros
X   = x - repmat(xmean,[T 1]);
Hx  = getHankelMatrix_MV(X,Ndelay);

Y = poolData(Hx(:,1:end-1)',Nvar,polyorder,usesine)'; % Take nonlinear measurement functions
Y = Y(2:end,:); % remove constant term
Yp = poolData(Hx(:,2:end)',n,polyorder,usesine)'; % Also for time-shifted data
Yp = Yp(2:end,:);   % remove constant term
[sysmodel_DMDc,U,Up] = DMDc(Y,Yp,Hu(1:end-1),dt);  % Use DMDc to estimate discrete-time, linear state-space model
sysmodel_DMDc.C(Nvar+1:end,Nvar+1:end) = 0;
Nt = length(t)-Ndelay+1;

%% Validation over training phase
Y0 = poolData(Hx(:,1)',Nvar,polyorder,usesine)';
[xDMDc,~] = lsim(sysmodel_DMDc,Hu,tspan(1:Nt),Y0(2:end));
xDMDc = xDMDc(:,1:Nvar);
xDMDc = xDMDc + repmat(xmean,[Nt 1]);

%% Show prediction over training data
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];
    ccolors(2,:)-[0.1 0.2 0.09];
    ccolors(3,:)-[0.1 0.2 0.09];
    ccolors(4,:)-[0.1 0.1 0.2];
    ccolors(5,:)-[0.1 0.2 0.09]];
for i = 1:Nvar
    ph(i) = semilogy(tspan,x(:,i),'-','Color',ccolors(i,:),'LineWidth',1); hold on
end
for i = 1:Nvar
    ph(Nvar+i) = semilogy(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
end
xlabel('Time')
ylabel('xi')
legend(ph([1,6]),'True',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);


%% Save Data
Model.name = ModelName;
Model.sys = sysmodel_DMDc;
Model.Ndelay = Ndelay;
Model.dt = dt;
Model.xmean = xmean;
Model.polyorder = polyorder;
Model.usesine = usesine;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')