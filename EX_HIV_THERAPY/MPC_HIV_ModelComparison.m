% LORENZ system
%118.2675
%98.2860
%2.1048e+03

clear all, close all, clc

SystemModel = 'HIV';

figpath = ['../../FIGURES/',SystemModel,'/'];mkdir(figpath)
datapath = ['../../DATA/',SystemModel,'/'];mkdir(datapath)
% datapath = ['/Users/ekaiser/Documents/Academia/Papers/KaKuBr_SINDYc-MPC/DATA/'];
addpath('../utils');

%% Load Models
Nvar = 5;
% InputSignalTypeModel = 'sine2'; % prbs; chirp; noise; sine2; sphs; mixed
% InputSignalTypeModel_Validation = 'sine3';

InputSignalTypeModel = 'prbs'; % prbs; chirp; noise; sine2; sphs; mixed
InputSignalTypeModel_Validation = 'prbs';%'sine2';

% DelayDMDc
ModelTypeDMDc = 'DMDc'; % actual
load(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelTypeDMDc,'_',InputSignalTypeModel,'.mat'])) % Model: DelayDMDc , DMDc
Models.DelayDMDc = Model;


% SINDYc
load(fullfile(datapath,['EX_',SystemModel,'_SI_SINDYc','_',InputSignalTypeModel,'.mat'])) % Model
Models.SINDYc = Model;

% NARX
% InputSignalTypeModel_NARX = 'sphs';
InputSignalTypeModel_NARX = InputSignalTypeModel;
% NARX_SUBSTRACT_MEAN = 1;
load(fullfile(datapath,['EX_',SystemModel,'_SI_NARX','_',InputSignalTypeModel_NARX,'.mat'])) % Model
Models.NARX = Model;

%% ADDITIONAL MODELS
ModelTypeDMDc = 'DelayDMDc';
load(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelTypeDMDc,'_',InputSignalTypeModel,'.mat'])) % Model: DelayDMDc , DMDc
Models.DelayDMDc2 = Model;

load(fullfile(datapath,['EX_',SystemModel,'_SI_partialSINDYc','_',InputSignalTypeModel,'_123.mat'])) % Model
Models.partialSINDYc = Model;

load(fullfile(datapath,['EX_',SystemModel,'_SI_eDMDc','_',InputSignalTypeModel,'.mat'])) % Model
Models.eDMDc = Model;
%% Get training data
ONLY_TRAINING_LENGTH = 1;
InputSignalType = InputSignalTypeModel;% prbs; chirp; noise; sine2; sphs; mixed
% Ndelay = Models.DelayDMDc.Ndelay;
getTrainingData
u = u';
%% Get validation data
InputSignalType = InputSignalTypeModel_Validation; %'other';% prbs; chirp; noise; sine2; sphs; mixed
% Ndelay = Models.DelayDMDc.Ndelay;
% getValidationData

%% FIGURE 1:  // Model comparison + Validation
tspanV =[tspan(end):dt:300];

% Forcing
if strcmp(InputSignalTypeModel_Validation,'sine2')
    forcing = @(x,t) forcing_sine2_positive(x,t);
    u_valid = forcing(0,tspanV')';
elseif strcmp(InputSignalTypeModel_Validation,'prbs')
    A = 1;
    taulim = [0.2 8];
    states = [0,0:0.25:1,0,0,0];
    Nswitch = 200;
    forcing = @(x,t) [A(1)*prbs(taulim, Nswitch, states, t,0,14)];
    u_valid = zeros(1,length(tspanV));
    for i = 1:length(tspanV)
        u_valid(:,i) = forcing(0,tspanV(i));
    end
    figure,plot(u_valid)
end

% to = t; clear t
% Reference
% x0 = x(1,1:Nvar);
x0 = [10, 0.1, 0.1, 1, 0.1];
[tA,xA] = ode45(@(t,x)HIVsys_ZURAKOWSKI(t,x,forcing(x,t)),tspanV,x0,options);   % true model

% DelayDMDc / DMDc
Ndelay = Models.DelayDMDc.Ndelay;
if Ndelay == 1
    %     x0      = [x(end,1:Nvar)];
    Hunew   = [u(end),u_valid(1:end-1)];
    [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspanV,[x0-[Models.DelayDMDc.xmean]]');
elseif Ndelay > 1
    %     x0      = [x(end-Ndelay+1,1:Nvar),x(end,1:Nvar)];
    Hunew   = [ u(end-Ndelay+1:end),u_valid(1:end-Ndelay);
        u(end),u_valid(1:end-1)];
    [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspanV,[x0-[Models.DelayDMDc.xmean]]');
    xB = xB(:,Nvar+1:2*Nvar); 
    %xB = xB + repmat(xmean,[size(xB,1) 1]);
end
xB = xB + repmat(Models.DelayDMDc.xmean,[length(tB) 1]);

% SINDYc
[tC,xC]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Models.SINDYc.Xi(:,1:Nvar),Models.SINDYc.polyorder,Models.SINDYc.usesine),tspanV,x0,options);  % approximate

% NARX
Udummy = [u(end),u_valid(1:end)];
Hdummy = zeros(Nvar,size(Udummy,2));
NARX_xmean = Models.NARX.xmean;
NARX_xnorm = Models.NARX.xnorm;

if Models.NARX.TRANSFORM_LOG
    x0_NARX = log(x0');
end
Hdummy(:,1) = x0_NARX-NARX_xmean;
Hdummy(:,1) = Hdummy(:,1)./NARX_xnorm;
[Us,Ui,Si] = preparets(Models.NARX.net,con2seq(Udummy),{},con2seq(Hdummy));
xD = Models.NARX.net(Us,Ui,Si); % Predict on validation data
xD = cell2mat(xD)';
if Models.NARX.TRANSFORM_LOG
    xD = exp(xD);
end
xD = [x0;xD(1:end-1,:)];
xD(2:end,:) = xD(2:end,:) + repmat(NARX_xmean',[size(xD,1)-1 1]);
xD(2:end,:) = xD(2:end,:).*repmat(NARX_xnorm',[size(xD,1)-1 1]);

%%
% partialSINDYc
options = odeset('RelTol',1e-14,'AbsTol',1e-14*ones(1,length(Models.partialSINDYc.SelectVars)));
[tE,xE]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Models.partialSINDYc.Xi(:,1:end-1),Models.partialSINDYc.polyorder,Models.partialSINDYc.usesine),tspanV,x0(Models.partialSINDYc.SelectVars),options);  % approximate
options = odeset('RelTol',1e-14,'AbsTol',1e-14*ones(1,n));

%%
% DelayDMDc / DMDc
Ndelay = Models.DelayDMDc2.Ndelay;
if Ndelay == 1
    Hunew   = [u(end),u_valid(1:end-1)];
    [xF,tF] = lsim(Models.DelayDMDc2.sys,Hunew,tspanV,[x0-[Models.DelayDMDc2.xmean]]');
elseif Ndelay > 1
    Hunew   = getHankelMatrix_MV(u_valid',Ndelay)';
    H0  = getHankelMatrix_MV(xA(1:Ndelay,:)-repmat(Models.DelayDMDc2.xmean,[Ndelay 1]),Ndelay);
    [xF,tF] = lsim(Models.DelayDMDc2.sys,Hunew,tspanV(Ndelay:end),H0);
    xF = xF(:,end-Nvar+1:end);
end
xF = xF + repmat(Models.DelayDMDc2.xmean,[length(tF) 1]);
xF = [nan(size(xA(1:Ndelay-1,:)));xF]; tF = [tA(1:Ndelay-1);tF];
%% eDMDc
Y0 = poolData([x0-Models.eDMDc.xmean],Nvar,Models.eDMDc.polyorder,Models.eDMDc.usesine)';
Hunew   = [u(end),u_valid(1:end-1)];

[xG,tG] = lsim(Models.eDMDc.sys,Hunew,tspanV,Y0(2:end));
xG = xG(:,1:Nvar);
xG = xG + repmat(Models.eDMDc.xmean,[length(tG) 1]);


%% Show results
clear ph
h = figure;
subplot(3,2,1), box on,
plot([tB(1) tB(1)],[1e-10 10],'--k'); hold on
plot(t,[x(:,1)],'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,1),'k','LineWidth',1.5);
plot(tB,xB(:,1),'r-','LineWidth',1.5);
plot(tC,xC(:,1),'g--','LineWidth',1.5);
plot(tA(1:end),xD(:,1),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
plot(tE,xE(:,1),'-','Color',[0,0,1],'LineWidth',1.5);
plot(tF,xF(:,1),':','Color',[0.5,0.3,.1],'LineWidth',1.5);
plot(tG,xG(:,1),'-.','Color',[1,0.3,0],'LineWidth',1.5);
grid on
axis tight
ylabel('x_1','FontSize',13)
set(gca,'FontSize',13)

subplot(3,2,2), box on
plot([tB(1) tB(1)],[1e-10 10],'--k'), hold on
plot(t,x(:,2),'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,2),'k','LineWidth',1.5);
plot(tB,xB(:,2),'r-','LineWidth',1.5);
plot(tC,xC(:,2),'g--','LineWidth',1.5);
plot(tA(1:end),xD(:,2),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
plot(tE,xE(:,2),'-','Color',[0,0,1],'LineWidth',1.5);
plot(tF,xF(:,2),':','Color',[0.5,0.3,0.1],'LineWidth',1.5);
plot(tG,xG(:,2),'-.','Color',[1,0.3,0],'LineWidth',1.5);
set(gca,'FontSize',13)
ylabel('x_2','FontSize',13)
axis tight

subplot(3,2,3), box on
semilogy([tB(1) tB(1)],[1e-1 100],'--k'), hold on
ph(1) = plot(t,x(:,3),'Color',[.4 .4 .4],'LineWidth',1.5);
ph(2) = plot(tA,xA(:,3),'k','LineWidth',1.5);
ph(3) = plot(tB,xB(:,3),'r-','LineWidth',1.5);
ph(4) = plot(tC,xC(:,3),'g--','LineWidth',1.5);
ph(5) = plot(tA(1:end),xD(:,3),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
plot(tE,xE(:,3),'-','Color',[0,0,1],'LineWidth',1.5);
plot(tF,xF(:,3),':','Color',[0.5,0.3,0.1],'LineWidth',1.5);
plot(tG,xG(:,3),'-.','Color',[1,0.3,0],'LineWidth',1.5);
grid on,axis tight
ylabel('x_3','FontSize',13)
set(gca,'FontSize',13)
xlabel('Time','FontSize',13)

% subplot(3,2,4), box on
% plot([tB(1) tB(1)],[1e-10 1],'--k'), hold on
% ph(1) = plot(t,x(:,4),'Color',[.4 .4 .4],'LineWidth',1.5);
% ph(2) = plot(tA,xA(:,4),'k','LineWidth',1.5);
% ph(3) = plot(tB,xB(:,4),'r-','LineWidth',1.5);
% ph(4) = plot(tC,xC(:,4),'g--','LineWidth',1.5);
% ph(5) = plot(tA(1:end),xD(:,4),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
% grid on,axis tight
% ylabel('x_3','FontSize',13)
% set(gca,'FontSize',13)
% xlabel('Time','FontSize',13)
% 
% subplot(3,2,5), box on
% plot([tB(1) tB(1)],[1e-10 8],'--k'), hold on
% ph(1) = plot(t,x(:,5),'Color',[.4 .4 .4],'LineWidth',1.5);
% ph(2) = plot(tA,xA(:,5),'k','LineWidth',1.5);
% ph(3) = plot(tB,xB(:,5),'r-','LineWidth',1.5);
% ph(4) = plot(tC,xC(:,5),'g--','LineWidth',1.5);
% ph(5) = plot(tA(1:end),xD(:,5),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
% l1=legend(ph,'Training','Validation',ModelTypeDMDc,'SINDYc','NARX');
% set(l1,'Location','NorthEast')
% l1.Position = [0.6 0.125 0.2 0.2];
% grid on,axis tight
% ylabel('x_3','FontSize',13)
% set(gca,'FontSize',13)
% xlabel('Time','FontSize',13)

set(h,'Units','Inches');
set(gcf,'Position',[1 1 6. 5.5])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
print(h,'-painters','-depsc2',  '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'.eps'],'-r0');

%% Show results
% tplot = (tA-tA(1))/7;
% for i = 1:3
% clear ph
% h = figure; box on, hold on, grid off
% plot(tplot,xA(:,i),'','Color',[0 0 0],'LineWidth',2);
% plot(tplot,xB(:,i),'r-','LineWidth',2);
% plot(tplot,xC(:,i),'g--','LineWidth',2);
% plot(tplot,xD(:,i),'-.','Color',[0.7,0.7,1],'LineWidth',2);
% axis tight
% xlabel('Time [weeks]'), ylabel(['x',num2str(i)])
% set(gca,'LineWidth',1, 'FontSize',16)
% set(gcf,'Position',[100 100 400 110])
% set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_x',num2str(i),'.eps']);
% end

ccols = [0.9,0.5,.5; 1,0,0; 0.5,0.15,0.15; 0.3,0.5,0.3; 0,1,0; 0.7,0.7,1];
lstyles = {'-','-','-','--','--','-.'};
new_order = [2,3,1,5,4,6];
tplot = (tA-tA(1))/7;
ylims = [-2,14;
         -3,6;
         -15,45;
         0,1;
         -2,5];
for i = 1:5
clear ph
h = figure; box on, hold on, grid off
plot(tplot,xA(:,i),'-','Color',[0 0 0],'LineWidth',1.5);
plot(tplot,xB(:,i),'Color',ccols(2,:),'LineWidth',1.5,'LineStyle',lstyles{2});
plot(tplot,xC(:,i),'Color',ccols(5,:),'LineWidth',1.5,'LineStyle',lstyles{5});
plot(tplot,xD(:,i),'Color',ccols(6,:),'LineWidth',1.5,'LineStyle',lstyles{6});
if i<4
    plot(tplot,xE(:,i),'Color',ccols(4,:),'LineWidth',1.5,'LineStyle',lstyles{4});
end
plot(tplot,xF(:,i),'Color',ccols(3,:),'LineWidth',1.5,'LineStyle',lstyles{3});
plot(tplot,xG(:,i),'Color',ccols(1,:),'LineWidth',1.5,'LineStyle',lstyles{1});
axis tight
ylim(ylims(i,:))
xlabel('Time [weeks]'), ylabel(['x',num2str(i)])
set(gca,'LineWidth',1, 'FontSize',16)
set(gcf,'Position',[100 100 400 130]) %110
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_x',num2str(i),'_COMP-ALL.eps']);
end

%% Error
errValid = zeros(Nmodels+1,length(tplot));

errValid(7,:) = sum(1/3*((xA(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2);
% errValid(2,:) = sum(1/3*(xB(:,1:3)-xA(:,1:3)).^2./repmat(max(xA(:,1:3),[],1),[length(tplot) 1]),2); %DMDc
% errValid(5,:) = sum(1/3*(xC(:,1:3)-xA(:,1:3)).^2./repmat(max(xA(:,1:3),[],1),[length(tplot) 1]),2); %SINDYc
% errValid(6,:) = sum(1/3*(xD(:,1:3)-xA(:,1:3)).^2./repmat(max(xA(:,1:3),[],1),[length(tplot) 1]),2); %NN
% errValid(4,:) = sum(1/3*(xE(:,1:3)-xA(:,1:3)).^2./repmat(max(xA(:,1:3),[],1),[length(tplot) 1]),2); %partialDSINDYc
% errValid(3,:) = sum(1/3*(xF(:,1:3)-xA(:,1:3)).^2./repmat(max(xA(:,1:3),[],1),[length(tplot) 1]),2); % Delay DMDc
% errValid(1,:) = sum(1/3*(xG(:,1:3)-xA(:,1:3)).^2./repmat(max(xA(:,1:3),[],1),[length(tplot) 1]),2); %eDMDc
errValid(2,:) = sum(1/3*((xB(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %DMDc
errValid(5,:) = sum(1/3*((xC(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %SINDYc
errValid(6,:) = sum(1/3*((xD(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %NN
errValid(4,:) = sum(1/3*((xE(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %partialDSINDYc
errValid(3,:) = sum(1/3*((xF(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); % Delay DMDc
errValid(1,:) = sum(1/3*((xG(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %eDMDc


clear ph
h = figure; box on, grid off
for i = 1:size(errValid,1)-1
% plot(tplot,xA(:,i),'-','Color',[0 0 0],'LineWidth',2);
semilogy(tplot,errValid(new_order(i),:),'Color',ccols(new_order(i),:),'LineWidth',1.5,'LineStyle',lstyles{new_order(i)}); hold on
end
axis tight
ylim([1e-15 1e5])
xlabel('Time [weeks]'), ylabel(['Error'])
set(gca,'LineWidth',1, 'FontSize',16,'ytick',[1e-15,1e-10,1e-5,1e0,1e5],'yticklabel',{'-15','-10','-05','-00','+05'})
set(gcf,'Position',[100 100 400 140]) %110
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_ValidationError_',InputSignalTypeModel,'_',InputSignalType,'_x',num2str(i),'_COMP-ALL.eps']);

%%
clear ph
h = figure; box on, hold on, grid off
plot(tplot,u_valid,'','Color',[0 0 0],'LineWidth',1.5);
axis tight
xlabel('Time [weeks]'), ylabel('u')
set(gca,'LineWidth',1, 'FontSize',16)
set(gcf,'Position',[100 100 400 100])%80
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_u.eps']);

%% Actuation signal
t_valid = tspanV;
clear ph
figure,box on, hold on,
ccolors = get(gca,'colororder');
plot([tA(1),tA(1)],[0 1],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
plot(t,u,'-k','LineWidth',1);
plot(t_valid,u_valid,'-k','LineWidth',1);
grid off
% ylim([min([u u_valid])+0.05*min([u u_valid]) max([u u_valid])+0.05*max([u u_valid])])
xlim([0 t_valid(end)])
xlabel('Time')
ylabel('Input')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_Actuation.eps']);


%% Truth
run_HIV_params
yout = poolDataLIST({'x1','x2','x3','x4','x5','u'},Models.SINDYc.Xi,6,3,0);

% ZURAKOWSKI
xi_truth.term{1} = {'1', 'x1', 'x1x2', 'x1x2u'};
xi_truth.coeff{1} = [lambda1, -d, -alpha1, eta*alpha1];
xi_truth.term{2} = {'x2', 'x1x2', 'x2x4', 'x2x5', 'x1x2u'};
xi_truth.coeff{2} = [-a, alpha1, -p1, -p2, -eta*alpha1];
xi_truth.term{3} = {'x3', 'x2x3', 'x1x2x3'};
xi_truth.coeff{3} = [-b2, -c2*q, c2];
xi_truth.term{4} = {'x4', 'x2x4'};
xi_truth.coeff{4} = [-b1, c1];
xi_truth.term{5} = {'x5', 'x2x3'};
xi_truth.coeff{5} = [-h, c2*q];

% Construct true Xi
Xi0 = zeros(size(Models.SINDYc.Xi));
for i = 1:Nvar
    for j = 1:length(xi_truth.term{i})
        idx = find(strcmp(yout,xi_truth.term{i}(j)));
        Xi0(idx,i) = xi_truth.coeff{i}(j);
    end
end


%% MPC
% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
Nweeks = 100;
Duration = Nweeks*7;                 % Run for 'Duration' time units
Ton = tA(end)+1;                         % Time units when control is turned on
x0n = [10, 0.1, 0.1, 0.1, 0.1]'; % Initial condition
Ts  = 1/24;             % Sampling time [1/hr]
Tcontrol = tA(end);
getMPCparams

% Reference state, which shall be achieved
xref = zeros(1,5);
xref(2) = ((c2*(lambda1-d*q)-b2*alpha1) - sqrt((c2*(lambda1-d*q)-b2*alpha1)^2 - 4*alpha1*c2*q*d*b2))/(2*alpha1*c2*q);
xref(1) = lambda1/(d+alpha1*xref(2));
xref(4) = 0;
xref(5) = (xref(2)*c2*(alpha1*q-a) + b2*alpha1)/(c2*p2*xref(2));
xref(3) = h*xref(5)/(c2*q*xref(2));


% Parameters Models
% ModelCollection = {ModelTypeDMDc, 'SINDYc', 'NARX'};
ModelCollection = {'eDMDc','DMDc', 'DelayDMDc', 'partialSINDYc', 'SINDYc', 'NARX'};
Nmodels = length(ModelCollection);

clear Results
Results(1:Nmodels) = struct('x',[], 'u', [], 't', [], 'xref', [], 'J', [], 'elapsed_time', []);

Ts = 1/12;%Models.SINDYc.dt; %dt;

%%
clear pest
for iM = Nmodels
    select_model = ModelCollection{iM}; % DelayDMDc; DMDc; NARX ; SINDYc
    
    % Prepare variables
    Nt = (Duration/Ts)+1;
    uopt0    = 0;
    xhat     = x0n;
    uopt     = uopt0.*ones(Nu,1);
    xHistory = zeros(Nvar,Nt); xHistory(:,1) = xhat;
    uHistory = zeros(1,Nt); uHistory(1)   = uopt(1);
    tHistory = zeros(1,Nt); tHistory(1)   = Tcontrol;
    rHistory = zeros(Nvar,Nt);
    rHistory(:,1) = xref;
    
    % Parameters Model
    switch select_model
        case 'eDMDc'
            pest.dt = Models.eDMDc.dt;
            pest.sys = Models.eDMDc.sys;
            pest.xmean = Models.eDMDc.xmean';
            pest.Nxlim = size(x,1);
            pest.polyorder = Models.eDMDc.polyorder;
            pest.usesine = Models.eDMDc.usesine;
        case 'DelayDMDc'
            pest.dt     = Models.DelayDMDc2.dt;
            pest.sys    = Models.DelayDMDc2.sys;
            pest.xmean  = Models.DelayDMDc2.xmean';
            pest.udelay = zeros(1,Models.DelayDMDc2.Ndelay);
%             pest.xdelay = zeros(Nvar,Ndelay);%[x(end-Ndelay+1,1:2)]';
            pest.xdelay = zeros(Nvar*Models.DelayDMDc2.Ndelay,1); % x0
            pest.Nxlim  = size(x,1);
            pest.Ndelay = Models.DelayDMDc2.Ndelay;
        case 'DMDc'
            pest.dt = Models.DelayDMDc.dt;
            pest.sys = Models.DelayDMDc.sys;
            pest.xmean = Models.DelayDMDc.xmean';
            pest.Nxlim = size(x,1);
        case 'SINDYc'
            pest.ahat = Models.SINDYc.Xi(:,1:Nvar); %Xi0(:,1:Nvar);%
            pest.polyorder = Models.SINDYc.polyorder;
            pest.usesine = Models.SINDYc.usesine;
            pest.dt = Models.SINDYc.dt;
            pest.SelectVars = 1:Nvar;
        case 'partialSINDYc'
            pest.ahat = Models.partialSINDYc.Xi(:,1:end-1); 
            pest.polyorder = Models.partialSINDYc.polyorder;
            pest.usesine = Models.partialSINDYc.usesine;
            pest.dt = Models.partialSINDYc.dt;  
            pest.SelectVars = Models.partialSINDYc.SelectVars;
        case 'NARX'
            pest.net = Models.NARX.net;
            pest.Ndelay = 1;
            pest.udelay = zeros(1,Ndelay);
            pest.xdelay = zeros(2,length(pest.udelay));
            pest.xdelay(:,1:Ndelay) = [x(end-Ndelay+1:end,1:2)]';
            pest.Nxlim = size(x,1);
            pest.TRANSFORM_LOG = Models.NARX.TRANSFORM_LOG;
    end
    
    
    % Start simulation
    fprintf('Simulation started.  It might take a while...\n')
    tic
    for ct = 1:(Duration/Ts)
        if (ct*Ts+Tcontrol)>Ton
        if mod(ct*Ts,7) == 0
            % NMPC with full-state feedback
            COSTFUN = @(u) ObjectiveFCN_models(u,xhat,N,Nu,xref,uHistory(:,ct),pest,diag(Q),R,Ru,select_model);
            CONSFUN = @(u) ConstraintFCN_models(u,uHistory(:,ct),xhat,N,LBo,UBo,LBdu,UBdu,pest,select_model);
            uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
        end
        end
        % Integrate system
        xhat = rk4u(@HIVsys_ZURAKOWSKI,xhat,uopt(1),Ts/2,2,[],0); %10, 7000: weekly
        xHistory(:,ct+1) = xhat;
        uHistory(:,ct+1) = uopt(1);
        tHistory(:,ct+1) = ct*Ts+Tcontrol;
        rHistory(:,ct+1) = xref;
        
        if strcmp(select_model,'DelayDMDc') && ((ct*Ts+Tcontrol)>Ton-Ts)
            pest.dt     = Models.DelayDMDc2.dt;
            pest.sys    = Models.DelayDMDc2.sys;
            pest.xmean  = Models.DelayDMDc2.xmean';
            pest.udelay = uHistory(:,(ct+1)-Models.DelayDMDc2.Ndelay+1:(ct+1));
            H0 = xHistory(:,(ct+1)-Models.DelayDMDc2.Ndelay+1:(ct+1))';%[x(end-Ndelay+1,:)]';
            pest.xdelay = getHankelMatrix_MV(H0-repmat(Models.DelayDMDc2.xmean,[pest.Ndelay 1]),pest.Ndelay);
%             pest.Nxlim  = size(x,1);    
        end
        
        if mod(ct,1000) == 0
            disp(['PROGRESS: ',num2str(100*ct/(Duration/Ts)),'%'])
        end
    end
    
    fprintf('Simulation finished!\n')
    %
    tElapsed = toc
    
    Results(iM).eval_time = tElapsed;
    Results(iM).xref = rHistory;
    Results(iM).x = xHistory;
    Results(iM).u = uHistory;
    Results(iM).t = tHistory;
    Results(iM).J = evalObjectiveFCN_HIV(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);
    Results(iM).elapsed_time = tElapsed;
   
    VIZ_SI_Validation_MPC
end

%N=10: 120.9128, 96.8967, 2.1755e+03
%%



% if iM==Nmodels
%     VIZ_SI_Validation_MPC_ensemble(ct,t_valid,u_valid,xHistory,tHistory,uHistory,Results,t,tA,xA,xB,xC,xD,select_model,SystemModel,figpath,N,InputSignalTypeModel,Nmodels,ModelCollection)
% end
%% Unforced
[tU,xU] = ode45(@(t,x)HIVsys_ZURAKOWSKI(t,x,0,[]),tHistory,x0n,options);   % true model
%%
xnormfinal = max(Results(5).x(:,1:ct),[],2);
tnorm = 7;
%%
for iM = 1:Nmodels
clear ph
figure;box on,hold on
ccolors = get(gca,'colororder');
plot(tHistory(1:ct)/tnorm,xref(1)*ones(length(tHistory(1:ct)),1)./xnormfinal(1),'--','Color',ccolors(1,:),'LineWidth',1), hold on
plot(tHistory(1:ct)/tnorm,xref(2)*ones(length(tHistory(1:ct)),1)./xnormfinal(2),'--','Color',ccolors(2,:),'LineWidth',1)
plot(tHistory(1:ct)/tnorm,xref(3)*ones(length(tHistory(1:ct)),1)./xnormfinal(3),'--','Color',ccolors(3,:),'LineWidth',1)
ph(1) = plot(Results(iM).t(1:ct)/tnorm,Results(iM).x(1,1:ct)./xnormfinal(1),'-','Color',ccolors(1,:),'LineWidth',1.5);
ph(2) = plot(Results(iM).t(1:ct)/tnorm,Results(iM).x(2,1:ct)./xnormfinal(2),'-','Color',ccolors(2,:),'LineWidth',1.5);
ph(3) = plot(Results(iM).t(1:ct)/tnorm,Results(iM).x(3,1:ct)./xnormfinal(3),'-','Color',ccolors(3,:),'LineWidth',1.5);
ph(4) = plot(Results(iM).t(1:ct)/tnorm,Results(iM).u(1:ct),'-k','LineWidth',1.5);
ylim([0 1])
legend(ph,'x1','x2','x3','u')
end

clear ph
figure;box on,hold on
ccolors = get(gca,'colororder');
plot(tHistory(1:ct)/tnorm,xref(1)*ones(length(tHistory(1:ct)),1)./xnormfinal(1),'--','Color',ccolors(1,:),'LineWidth',1), hold on
plot(tHistory(1:ct)/tnorm,xref(2)*ones(length(tHistory(1:ct)),1)./xnormfinal(2),'--','Color',ccolors(2,:),'LineWidth',1)
plot(tHistory(1:ct)/tnorm,xref(3)*ones(length(tHistory(1:ct)),1)./xnormfinal(3),'--','Color',ccolors(3,:),'LineWidth',1)
ph(1) = plot(tU/7,xU(:,1)./xnormfinal(1),'-','Color',ccolors(1,:),'LineWidth',1.5);
ph(2) = plot(tU/7,xU(:,2)./xnormfinal(2),'-','Color',ccolors(2,:),'LineWidth',1.5);
ph(3) = plot(tU/7,xU(:,3)./xnormfinal(3),'-','Color',ccolors(3,:),'LineWidth',1.5);
ylim([0 1])
legend(ph,'x1','x2','x3')

%%
% figure,hold on
% plot([1 length(Results(1).J); 1 length(Results(1).J); 1 length(Results(1).J)]',[xref1,xref1]','-k')
% plot(Results(1).x','-'), plot(Results(2).x','--')
%
% figure,plot(cumsum(Results(1).J'),'-'), hold on, plot(cumsum(Results(2).J'),'--')
% figure,plot(Results(1).u','-'), hold on, plot(Results(2).u','--')

%%

% iM = 1; [J1,Js1,Ju1] = evalObjectiveFCN(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);
%
% iM = 2; [J2,Js2,Ju2] = evalObjectiveFCN(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);
%
% figure,
% subplot(3,1,1), plot(cumsum(J1),'-k'), hold on, plot(cumsum(J2),'--r')
% subplot(3,1,2), plot(cumsum(Js1),'-k'), hold on, plot(cumsum(Js2),'--r')
% subplot(3,1,3), plot(cumsum(Ju1),'-k'), hold on, plot(cumsum(Ju2),'--r')
% return
% figure,plot(cumsum(Results(1).J'),'-'), hold on, plot(cumsum(Results(2).J'),'--')
% plot(cumsum(J'),'-.k')
%
% figure; hold on, grid on
% plot(Results(1).x'-Results(iM).xref','-k'), plot(Results(2).x'-Results(iM).xref','--r')
%
% figure,plot(cumsum(Results(1).J'),'-k'), hold on, plot(cumsum(Results(2).J'),'--r')
% plot(cumsum(Js'),'-.g')
%
% figure,plot((Results(1).J'),'-k'), hold on, plot((Results(2).J'),'--r')
% plot((Js'),'-.g')
% return
% iM = 2; [J,Js,Ju] = evalObjectiveFCN(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);



%% Show cost
% ModelCollection = {'eDMDc','DMDc', 'DelayDMDc', 'partialSINDYc', 'SINDYc', 'NARX'};
ccols = [0.9,0.5,.5; 1,0,0; 0.5,0.15,0.15; 0.3,0.5,0.3; 0,1,0; 0.7,0.7,1];
lstyles = {'-','-','-','--','--','-.'};
new_order = [2,3,1,5,4,6];
clear ph
figure, hold on, box on
% ph(2) = plot(Results(1).t./7,cumsum(Results(2).J),'-r','LineWidth',2);
% ph(5) = plot(Results(1).t./7,cumsum(Results(5).J),'--g','LineWidth',2);
% ph(6) = plot(Results(1).t./7,cumsum(Results(6).J),'-.','Color',[0.7,0.7,1],'LineWidth',2);
% 
% ph(1) = plot(Results(1).t./7,cumsum(Results(1).J),'-','Color',[0.9,0.5,.5],'LineWidth',2);
% ph(3) = plot(Results(1).t./7,cumsum(Results(3).J),'-','Color',[0.5,0.15,0.15],'LineWidth',2);
% ph(4) = plot(Results(1).t./7,cumsum(Results(4).J),'--','Color',[0.3,0.5,0.3],'LineWidth',2);

for iM = 1:Nmodels
    ph(iM) = plot(Results(1).t./7,cumsum(Results(new_order(iM)).J),'-','Color',ccols(new_order(iM),:), ...
        'LineWidth',2,'LineStyle',lstyles{new_order(iM)});
end

xlim([Results(1).t(1)./7 Results(1).t(end)./7])
% ylim([0 4.5])
ylabel('Cost','FontSize',14)
xlabel('Time [weeks]','FontSize',14)
% set(gca,'yscale','log')
set(gca,'LineWidth',1, 'FontSize',14)
% set(gcf,'Position',[100 100 300 200])
set(gcf,'Position',[100 100 300 300])
set(gcf,'PaperPositionMode','auto')

print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_COMP-ALL.eps']);

l1=legend(ph,ModelCollection(new_order));
set(l1,'Location','NorthWest','FontSize',20)
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_COMP-ALL_wLegend.eps']);


axis([0 1 max(cumsum(Results(new_order(iM)).J)) max(cumsum(Results(new_order(iM)).J))+10])
axis off
set(l1,'Location','SouthWest','FontSize',20)
l1.Position = [0.2 0.2 l1.Position(3:4)];
set(gcf,'Position',[100 100 200 200])
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_COMP-ALL_Legend.eps']);


%% Show x1-3 for each model
% ccols = [1 0 0; 0.7,0.7,1; 0 1 0];
ccols = zeros(Nmodels,3);

for ivar = Nmodels
    ylimits = [0 11; 0 9; 0 1500];
    clear ph hh
    figure, hold on, box on
    hh = area(Results(1).t./7,Results(ivar).u,'LineStyle','none');
    hh.FaceColor = 0.8*ones(1,3);
    ph(4) = plot(Results(1).t(2:end)./7,(Results(ivar).xref(1,2:end))./xnormfinal(1),'-k','LineWidth',1);
    ph(5) = plot(Results(1).t(2:end)./7,(Results(ivar).xref(2,2:end))./xnormfinal(2),':k','LineWidth',1);
    ph(6) = plot(Results(1).t(2:end)./7,(Results(ivar).xref(3,2:end))./xnormfinal(3),'-.k','LineWidth',1);
    ph(1) = plot(Results(1).t./7,(Results(ivar).x(1,:))./xnormfinal(1),'-','Color',ccols(ivar,:),'LineWidth',2);
    ph(2) = plot(Results(1).t./7,(Results(ivar).x(2,:))./xnormfinal(2),':','Color',ccols(ivar,:),'LineWidth',2);
    ph(3) = plot(Results(1).t./7,(Results(ivar).x(3,:))./xnormfinal(3),'-.','Color',ccols(ivar,:),'LineWidth',2);
    ylim([0 1.2])
    xlim([min(Results(1).t./7) max(Results(1).t./7)])
    ylabel('xi','FontSize',14)
    xlabel('Time [weeks]','FontSize',14)
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 300])
    %     set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_',ModelCollection{ivar},'_COMP-ALL.eps']);  
    
end

    axis off
    l1=legend(ph,'-x1--','-x2--','-x3--','-Ref--');
    set(l1,'Location','SouthEast','Orientation','horizontal','FontSize',14)
    l1.Position = [l1.Position(1)+0.1 l1.Position(2)+0.1 l1.Position(3:4)];
    xlim([Results(1).t(1) Results(1).t(end)])
    % ylim([10 1500])
    ylabel('xi','FontSize',14)
    xlabel('Time','FontSize',14)
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 50])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_Model_legend_COMP-ALL.eps']);

    %% No control
    
    ylimits = [0 11; 0 9; 0 1500];
    clear ph hh
    figure, hold on, box on
%     hh = area(Results(1).t./7,Results(ivar).u,'LineStyle','none');
%     hh.FaceColor = 0.8*ones(1,3);
    ph(4) = plot(Results(1).t(2:end)./7,(Results(ivar).xref(1,2:end))./xnormfinal(1),'-k','LineWidth',1);
    ph(5) = plot(Results(1).t(2:end)./7,(Results(ivar).xref(2,2:end))./xnormfinal(2),':k','LineWidth',1);
    ph(6) = plot(Results(1).t(2:end)./7,(Results(ivar).xref(3,2:end))./xnormfinal(3),'-.k','LineWidth',1);
    ph(1) = plot(Results(1).t./7,(xU(:,1))./xnormfinal(1),'-','Color','k','LineWidth',2);
    ph(2) = plot(Results(1).t./7,(xU(:,2))./xnormfinal(2),':','Color','k','LineWidth',2);
    ph(3) = plot(Results(1).t./7,(xU(:,3))./xnormfinal(3),'-.','Color','k','LineWidth',2);
    ylim([0 1.2])
    xlim([min(Results(1).t./7) max(Results(1).t./7)])
    ylabel('xi','FontSize',14)
    xlabel('Time [weeks]','FontSize',14)
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 300])
    %     set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_UNFORCED','_COMP-ALL.eps']);  
  
    %%
    set(gcf,'Position',[100 100 300 200])
    legend off
    axis on
    xlim([Results(1).t(end) Results(1).t(end)+1])
    ylim([10000 20000])
    l1=legend(ph,'-x1--','-x2--','-x3--','-Ref--');
    set(l1,'Location','SouthEast','Orientation','vertical','FontSize',14)
%     l1.Position = [l1.Position(1)+0.1 l1.Position(2)+0.1 l1.Position(3:4)];
    axis off
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 120 150])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_Model_legend_COMP-ALL_short.eps']);

    %%
    set(gcf,'Position',[100 100 300 200])
    legend off
    axis on
    xlim([Results(1).t(end) Results(1).t(end)+1])
    ylim([10000 20000])
    l1=legend(ph,'-x1--','-x2--','-x3--','-Ref1--','-Ref2--','-Ref3--');
    set(l1,'Location','SouthEast','Orientation','vertical','FontSize',14)
%     l1.Position = [l1.Position(1)+0.1 l1.Position(2)+0.1 l1.Position(3:4)];
    axis off
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 120 160])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_Model_legend_COMP-ALL_long.eps']);
    
%% Show control stage separate
% clear ph
% figure, hold on, box on
% ph(4) = plot(Results(1).t(2:end),(Results(1).xref(1,2:end)),'-k','LineWidth',2);
% ph(1) = plot(Results(1).t,(Results(1).x(1,:)),'-r','LineWidth',2);
% ph(3) = plot(Results(1).t,(Results(3).x(1,:)),'-.','Color',[0.7,0.7,1],'LineWidth',2);
% ph(2) = plot(Results(1).t,(Results(2).x(1,:)),'--g','LineWidth',2);
%
% l1=legend(ph,ModelCollection,'Ref');
% set(l1,'Location','NorthEast')
% xlim([Results(1).t(1) Results(1).t(end)])
% % ylim([-0.25 .25])
% ylabel('x1','FontSize',14)
% xlabel('Time','FontSize',14)
% % set(gca,'yscale','log')
% set(gca,'LineWidth',1, 'FontSize',14)
% set(gcf,'Position',[100 100 300 200])
% set(gcf,'PaperPositionMode','auto')

return
% print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_x1','_N_',num2str(N),'.eps']);


%%
clear ph
figure, hold on, box on
ph(1) = plot(Results(1).t,(Results(1).u(:)),'-r','LineWidth',2);
ph(3) = plot(Results(1).t,(Results(3).u(:)),'-.','Color',[0.7,0.7,1],'LineWidth',2);
ph(2) = plot(Results(1).t,(Results(2).u(:)),'--g','LineWidth',2);

l1=legend(ph,ModelCollection,'Ref');
set(l1,'Location','NorthEast')
xlim([Results(1).t(1) Results(1).t(end)])
% ylim([-0.25 .25])
ylabel('Input','FontSize',14)
xlabel('Time','FontSize',14)
% set(gca,'yscale','log')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

% print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Input','_N_',num2str(N),'.eps']);
%%
z=1

%% Save results
save(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_N_',num2str(N),'_COMP-ALL.mat']),'Results')
save(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_All','_N_',num2str(N),'_COMP-ALL.mat']))