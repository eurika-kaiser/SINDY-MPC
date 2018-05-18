% Automatic flight control of F8
% System identification: SINDYc

clear all, close all, clc
figpath = '../../FIGURES/F8/'; mkdir(figpath)
datapath = '../../DATA/F8/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'F8';
Nvar = 3;

%% Generate Data %{'sine2', 'chirp','prbs', 'sphs'}
InputSignalType = 'sine3';%'sine2'; sine3, %prbs; noise; sine2; sphs; 
ONLY_TRAINING_LENGTH = 1;
ENSEMBLE_DATA = 0;
getTrainingData
u = u';

%% SINDYc
% Parameters
ModelName = 'SINDYc'; 
polyorder = 3;
usesine = 0;
% lambda_vec = [.01,.01,0.1]; % sphs // WORKED IN MPC
lambda_vec = [.0001,.01,0.01]; % sine2 // WORKED IN MPC
eps = 0;

%Add noise
%eps = 0.05*std(x(:,1));
%x = x + eps*randn(size(x));       

trainSINDYc

%% Compare with true model parameters
Xi0 = zeros(size(Xi));
Xi0([2,4,5,6,8,10,16,18,19,25,35],1) = [-0.877,1,-0.215,0.47,-0.088,-0.019,3.846,-1,0.28,0.47,0.63];
Xi0([4],2) = 1;
Xi0([2,4,5,6,16,19,25,35],3) = [-4.208,-0.396,-20.967,-0.47,-3.564,6.265,46,61.4];
Xi-Xi0

%% Prediction  over training phase
if any(strcmp(InputSignalType,{'sine2', 'chirp','prbs', 'sphs'})==1)
    [tSINDYc,xSINDYc]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xi(:,1:Nvar),polyorder,usesine),tspan,x0,options);  % approximate
else
    p.ahat = Xi(:,1:Nvar);
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

%% Show validation
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ccolors_valid = [ccolors(1,:)-[0 0.2 0.2];ccolors(2,:)-[0.1 0.2 0.09];ccolors(3,:)-[0.1 0.2 0.09]];
for i = 1:Nvar
    ph(i) = plot(tspan,x(:,i),'-','Color',ccolors(i,:),'LineWidth',1); hold on
end
for i = 1:Nvar
    ph(Nvar+i) = plot(tspan,xSINDYc(:,i),'--','Color',ccolors_valid(i,:),'LineWidth',2);
end
ylim([-0.8 0.9])
xlabel('Time')
ylabel('xi')
legend(ph([1,4]),'True',ModelName)
legend(ph([1,2,3]),'angle of attack','pitch angle', 'pitch rate')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);

%% Prediction
% Truth
tspanV   = [100:dt:200];
xA      = xv;
tA      = tv;

% Model
if any(strcmp(InputSignalType,{'sine2', 'chirp','prbs', 'sphs'})==1)
    [tB,xB]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xi(:,1:Nvar),polyorder,usesine),tspanV,x(end,:),options);  % approximate
    xB = xB(2:end,:);
    tB = tB(2:end); 
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
u = u';
VIZ_SI_Validation

%% Save Data
Model.name = 'SINDYc';
Model.polyorder = polyorder;
Model.usesine = usesine;
Model.Xi = Xi;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')