% LOTKA-VOLTERRA system
% System identification: DelayDMDc

clear all, close all, clc
figpath = '../../FIGURES/LORENZ_SLIDES/'; mkdir(figpath)
% datapath = '../../DATA/LORENZ/';
addpath('../utils');


SystemModel = 'LORENZ';
ModelName = 'SINDYc';
%% Generate Data %{'sine2', 'chirp','prbs', 'sphs'}
InputSignalType = 'sphs'; %prbs; chirp; noise; sine2; sphs; mixed
ONLY_TRAINING_LENGTH = 1;
getTrainingData

%% SINDYc
u = u';

% Parameters
ModelName = 'SINDYc'; iModel = 2;
Nvar = 3;
polyorder = 4;
usesine = 0;
lambda = 0.1;
eps = 0;
xclean = x;

%Add noise
xclean = x;
eps = 0.05*std(x(:,1));
x = x + eps*randn(size(x));

%%
figure;
plot3(x(:,1),x(:,2),x(:,3),'-k','LineWidth',1);
color_line3(xclean(:,1),xclean(:,2),xclean(:,3),u,'LineWidth',2);
colormap(parula)
view(30,15)
axis tight
axis off
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_PhasePlot.eps']);

%%
polyorder = 2;
usesine = 0;
lambda = 0.5;
trainSINDYc

return
%%
ii = 1;
figure,
plot(t(1:size(xaug(:,ii),1)),xaug(:,ii),'-r','LineWidth',2);
hold on
plot(t(1:size(xaug(:,ii),1)),xclean(50:end-50-1,ii),'-k','LineWidth',1)
ylim([min(xaug(:,ii))-0.5*xabs  max(xaug(:,ii))+0.5*xabs])
xlim([1 4])
xlabel('t'), ylabel(['x_',num2str(ii)])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 200 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_x',num2str(ii),'_eta',num2str(eps),'.eps']);

return
%%
ccol = {'b','r','g','k'};
tshort = t(1:size(dx(:,1),1));

for ii = 1:3

xabs = max(abs([min(xaug(:,ii)), max(xaug(:,ii))]));    
figure;
plot(tshort,xaug(:,ii),'-','Color',ccol{ii},'LineWidth',2)
ylim([min(xaug(:,ii))-0.5*xabs  max(xaug(:,ii))+0.5*xabs])
xlim([min(tshort)-0.05*max(tshort) max(t)+0.05*max(tshort)])
% daspect([1 30 1])
xlabel('t'), ylabel(['x_',num2str(ii)])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_x',num2str(ii),'_eta',num2str(eps),'.eps']);

axis off
axis tight
% daspect([1 50 1])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_x',num2str(ii),'_AxisOFF','_eta',num2str(eps),'.eps']);

% dx
xabs = max(abs([min(dx(:,ii)), max(dx(:,ii))]));  

figure;
plot(tshort,dx(:,ii),'-','Color',ccol{ii},'LineWidth',2)
ylim([min(dx(:,ii))-0.5*xabs  max(dx(:,ii))+0.5*xabs])
xlim([min(tshort)-0.05*max(tshort) max(t)+0.05*max(tshort)])
% daspect([1 30 1])
xlabel('t'), ylabel(['dx_',num2str(ii)])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_dx',num2str(ii),'_eta',num2str(eps),'.eps']);

axis off
axis tight
% daspect([1 50 1])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_dx',num2str(ii),'_AxisOFF','_eta',num2str(eps),'.eps']);
end

ii = 4; % u
xabs = max(abs([min(xaug(:,ii)), max(xaug(:,ii))]));    
figure;
plot(tshort,xaug(:,ii),'-','Color',ccol{ii},'LineWidth',2)
ylim([min(xaug(:,ii))-0.5*xabs  max(xaug(:,ii))+0.5*xabs])
xlim([min(tshort)-0.05*max(tshort) max(t)+0.05*max(tshort)])
% daspect([1 30 1])
xlabel('t'), ylabel(['u'])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_u','_eta',num2str(eps),'.eps']);

axis off
axis tight
% daspect([1 50 1])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_u','_AxisOFF','_eta',num2str(eps),'.eps']);


return
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
    xSINDYc = xSINDYc';%(:,2:N+1)';
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
xlim([0 (length(tspan)-1)*dt]), ylim([-25 50])
xlabel('Time')
ylabel('xi')
% legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
legend(ph([1,4]),'True',ModelName)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '-cmyk', [figpath,'EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.eps']);

%% Validation 3D
filename = ['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_train'];
xModel = xSINDYc;
xTRUTH = x;
color_type = 'models';
VIZ_3D_MODELvsTRUTH

%% Prediction
% Reference
tspanV   = [10:dt:20];
xA      = xv;
tA      = tv;

% Model
if any(strcmp(InputSignalType,{'sine2', 'chirp','prbs', 'sphs'})==1)
    [tB,xB]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xi(:,1:Nvar),polyorder,usesine),tspanV,x(end,:),options);  % approximate
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

%% Prediction
filename = ['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'_valid'];
xModel = xB;
xTRUTH = xA;
color_type = 'models';
VIZ_3D_MODELvsTRUTH

%% Show training and prediction
VIZ_SI_Validation

%% Error // TODO

%% Save Data
Model.name = 'SINDYc';
Model.polyorder = polyorder;
Model.usesine = usesine;
Model.Xi = Xi;
Model.dt = dt;
save(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelName,'_',InputSignalType,'.mat']),'Model')