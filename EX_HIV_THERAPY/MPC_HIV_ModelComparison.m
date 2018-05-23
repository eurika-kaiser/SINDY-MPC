% Optimal drug treatment of HIV using model predictive control
% and various models

clear all, close all, clc

SystemModel = 'HIV';

figpath = ['../FIGURES/',SystemModel,'/'];mkdir(figpath)
datapath = ['../DATA/',SystemModel,'/'];mkdir(datapath)
addpath('../utils');

%% Parameters
Nvar = 5;                                   % Number of variables
DATA_ENSEMBLE = 0;                          % Use data from single trajectory
InputSignalTypeModel = 'prbs';              % Actuation type used to identify mod
InputSignalTypeModel_Validation = 'prbs';
ModelCollection = {'eDMDc','DMDc', 'DelayDMDc', 'partialSINDYc', 'SINDYc', 'NARX'};
Nmodels = length(ModelCollection);
%% Load Models
% DelayDMDc
ModelTypeDMDc = 'DMDc';
load(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelTypeDMDc,'_',InputSignalTypeModel,'.mat'])) % Model: DelayDMDc , DMDc
Models.DelayDMDc = Model;

% SINDYc
load(fullfile(datapath,['EX_',SystemModel,'_SI_SINDYc','_',InputSignalTypeModel,'.mat'])) % Model
Models.SINDYc = Model;

% NARX
% NARX_SUBSTRACT_MEAN = 1;
load(fullfile(datapath,['EX_',SystemModel,'_SI_NARX','_',InputSignalTypeModel,'.mat'])) % Model
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
InputSignalType = InputSignalTypeModel;
getTrainingData
u = u';

%% Get validation data + results
InputSignalType = InputSignalTypeModel_Validation;
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

% Reference // True dynamics
x0 = [10, 0.1, 0.1, 1, 0.1];
[tA,xA] = ode45(@(t,x)HIVsys_ZURAKOWSKI(t,x,forcing(x,t)),tspanV,x0,options);   % true model

% DMDc
Hunew   = [0,u_valid(1:end-1)];
[xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspanV,[x0-[Models.DelayDMDc.xmean]]');
xB = xB + repmat(Models.DelayDMDc.xmean,[length(tB) 1]);

% SINDYc
[tC,xC]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Models.SINDYc.Xi(:,1:Nvar),Models.SINDYc.polyorder,Models.SINDYc.usesine),tspanV,x0,options);

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
xD(2:end,:) = xD(2:end,:) + repmat(NARX_xmean',[size(xD,1)-1 1]); % mean might be zero and has no effect

% partialSINDYc
options = odeset('RelTol',1e-14,'AbsTol',1e-14*ones(1,length(Models.partialSINDYc.SelectVars))); % different number of variables
[tE,xE]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Models.partialSINDYc.Xi(:,1:end-1),Models.partialSINDYc.polyorder,Models.partialSINDYc.usesine),tspanV,x0(Models.partialSINDYc.SelectVars),options);  % approximate
options = odeset('RelTol',1e-14,'AbsTol',1e-14*ones(1,n)); % change back

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

% eDMDc
Y0 = poolData([x0-Models.eDMDc.xmean],Nvar,Models.eDMDc.polyorder,Models.eDMDc.usesine)';
Hunew   = [u(end),u_valid(1:end-1)];

[xG,tG] = lsim(Models.eDMDc.sys,Hunew,tspanV,Y0(2:end));
xG = xG(:,1:Nvar);
xG = xG + repmat(Models.eDMDc.xmean,[length(tG) 1]);


%% Show validation results (only first 3 variables as partialSINDYc uses only those)
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

%% Show results // different plots

ccols = [0.9,0.5,.5; 1,0,0; 0.5,0.15,0.15; 0.3,0.5,0.3; 0,1,0; 0.7,0.7,1];
lstyles = {'-','-','-','--','--','-.'};
new_order = [2,3,1,5,4,6];
tplot = (tA-tA(1))/7; % shift and rescale time axis (weekly, starting at 0)
ylims = [-2,14;
    -3,6;
    -15,45;
    0,1;
    -2,5];
for i = 1:Nvar
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

%% Compute and plot Error
errValid = zeros(Nmodels+1,length(tplot));

errValid(7,:) = sum(1/3*((xA(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2);
errValid(2,:) = sum(1/3*((xB(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %DMDc
errValid(5,:) = sum(1/3*((xC(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %SINDYc
errValid(6,:) = sum(1/3*((xD(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %NN
errValid(4,:) = sum(1/3*((xE(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %partialDSINDYc
errValid(3,:) = sum(1/3*((xF(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); % Delay DMDc
errValid(1,:) = sum(1/3*((xG(:,1:3)-xA(:,1:3))./repmat(max(xA(:,1:3),[],1),[length(tplot) 1])).^2,2); %eDMDc


clear ph
h = figure; box on, grid off
for i = 1:size(errValid,1)-1
    semilogy(tplot,errValid(new_order(i),:),'Color',ccols(new_order(i),:),'LineWidth',1.5,'LineStyle',lstyles{new_order(i)}); hold on
end
axis tight
ylim([1e-15 1e5])
xlabel('Time [weeks]'), ylabel(['Error'])
set(gca,'LineWidth',1, 'FontSize',16,'ytick',[1e-15,1e-10,1e-5,1e0,1e5],'yticklabel',{'-15','-10','-05','-00','+05'})
set(gcf,'Position',[100 100 400 140]) %110
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_ValidationError_',InputSignalTypeModel,'_',InputSignalType,'_x',num2str(i),'_COMP-ALL.eps']);

%% Show applied actuation input
clear ph
h = figure; box on, hold on, grid off
plot(tplot,u_valid,'','Color',[0 0 0],'LineWidth',1.5);
axis tight
xlabel('Time [weeks]'), ylabel('u')
set(gca,'LineWidth',1, 'FontSize',16)
set(gcf,'Position',[100 100 400 100])%80
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_u.eps']);

%% Actuation signal over training and validation stage
t_valid = tspanV;
clear ph
figure,box on, hold on,
ccolors = get(gca,'colororder');
plot([tA(1),tA(1)],[0 1],':','Color',[0.4,0.4,0.4],'LineWidth',4)
plot(t,u,'-k','LineWidth',1);
plot(t_valid,u_valid,'-k','LineWidth',1);
grid off
xlim([0 t_valid(end)])
xlabel('Time')
ylabel('Input')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_u_train-valid.eps']);


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
Nweeks = 50;
DaysPerWeek  = 7;
Duration = Nweeks*DaysPerWeek;                 % Run for 'Duration' time units
Ton = tA(end)+1;                     % Time units when control is turned on
x0n = [10, 0.1, 0.1, 0.1, 0.1]';     % Initial condition
Ts = 1/12;                           % Sampling time [0.5/hr]
Tcontrol = tA(end);                  % Turn control on at time
getMPCparams

% Reference state, which shall be approached
xref = zeros(1,5);
xref(2) = ((c2*(lambda1-d*q)-b2*alpha1) - sqrt((c2*(lambda1-d*q)-b2*alpha1)^2 - 4*alpha1*c2*q*d*b2))/(2*alpha1*c2*q);
xref(1) = lambda1/(d+alpha1*xref(2));
xref(4) = 0;
xref(5) = (xref(2)*c2*(alpha1*q-a) + b2*alpha1)/(c2*p2*xref(2));
xref(3) = h*xref(5)/(c2*q*xref(2));


% Parameters Models: Select models for MPC
ModelCollection = {'eDMDc','DMDc', 'DelayDMDc', 'partialSINDYc', 'SINDYc', 'NARX'};
Nmodels = length(ModelCollection);

clear Results
Results(1:Nmodels) = struct('x',[], 'u', [], 't', [], 'xref', [], 'J', [], 'elapsed_time', []);

%% RUN MPC over selected models
clear pest
for iM = 1:1%Nmodels
    select_model = ModelCollection{iM};
    
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
            pest.xdelay = zeros(Nvar*Models.DelayDMDc2.Ndelay,1);
            pest.Nxlim  = size(x,1);
            pest.Ndelay = Models.DelayDMDc2.Ndelay;
        case 'DMDc'
            pest.dt = Models.DelayDMDc.dt;
            pest.sys = Models.DelayDMDc.sys;
            pest.xmean = Models.DelayDMDc.xmean';
            pest.Nxlim = size(x,1);
        case 'SINDYc'
            pest.ahat = Models.SINDYc.Xi(:,1:Nvar); %Xi0(:,1:Nvar); % Change to execute true model
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
    for ct = 1:(Duration/Ts) % Iterate over control stage
        if (ct*Ts+Tcontrol)>Ton  % Turn control on at a certain time
            if mod(ct*Ts,7) == 0  % Update (Optimize) control input once a week
                % Nonlinear MPC optimization
                COSTFUN = @(u) ObjectiveFCN_models(u,xhat,N,Nu,xref,uHistory(:,ct),pest,diag(Q),R,Ru,select_model);
                CONSFUN = @(u) ConstraintFCN_models(u,uHistory(:,ct),xhat,N,LBo,UBo,LBdu,UBdu,pest,select_model);
                uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
            end
        end
        
        % Integrate system
        xhat = rk4u(@HIVsys_ZURAKOWSKI,xhat,uopt(1),Ts/2,2,[],0); % use half the timestep and keep control constant for smaller error
        xHistory(:,ct+1) = xhat;
        uHistory(:,ct+1) = uopt(1);
        tHistory(:,ct+1) = ct*Ts+Tcontrol;
        rHistory(:,ct+1) = xref;
        
        if strcmp(select_model,'DelayDMDc') && ((ct*Ts+Tcontrol)>Ton-Ts) % Update history of states and inputs in model parameters
            pest.dt     = Models.DelayDMDc2.dt;
            pest.sys    = Models.DelayDMDc2.sys;
            pest.xmean  = Models.DelayDMDc2.xmean';
            pest.udelay = uHistory(:,(ct+1)-Models.DelayDMDc2.Ndelay+1:(ct+1));
            H0 = xHistory(:,(ct+1)-Models.DelayDMDc2.Ndelay+1:(ct+1))';
            pest.xdelay = getHankelMatrix_MV(H0-repmat(Models.DelayDMDc2.xmean,[pest.Ndelay 1]),pest.Ndelay);
        end
        
        if mod(ct,1000) == 0
            disp(['PROGRESS: ',num2str(100*ct/(Duration/Ts)),'%'])
        end
    end
    tElapsed = toc
    fprintf('Simulation finished!\n')
    
    
    % Collect results
    Results(iM).eval_time = tElapsed;
    Results(iM).xref = rHistory;
    Results(iM).x = xHistory;
    Results(iM).u = uHistory;
    Results(iM).t = tHistory;
    Results(iM).J = evalObjectiveFCN_HIV(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);
    Results(iM).elapsed_time = tElapsed;
    
end


%% Run unforced system
[tU,xU] = ode45(@(t,x)HIVsys_ZURAKOWSKI(t,x,0,[]),tHistory,x0n,options);   % true model

%% Normalize states and rescale time [days --> weeks]
xnormfinal = max(Results(1).x(:,1:ct),[],2);
tnorm = 7;

%% Show results
for iM = 1:1:Nmodels
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

%% Show cost
% ModelCollection = {'eDMDc','DMDc', 'DelayDMDc', 'partialSINDYc', 'SINDYc', 'NARX'};
ccols = [0.9,0.5,.5; 1,0,0; 0.5,0.15,0.15; 0.3,0.5,0.3; 0,1,0; 0.7,0.7,1];
lstyles = {'-','-','-','--','--','-.'};
new_order = [2,3,1,5,4,6]; % Show model results in different order
clear ph
figure, hold on, box on
for iM = 1:Nmodels
    ph(iM) = plot(Results(1).t./7,cumsum(Results(new_order(iM)).J),'-','Color',ccols(new_order(iM),:), ...
        'LineWidth',2,'LineStyle',lstyles{new_order(iM)});
end

xlim([Results(1).t(1)./7 Results(1).t(end)./7])
ylabel('Cost','FontSize',14)
xlabel('Time [weeks]','FontSize',14)
set(gca,'LineWidth',1, 'FontSize',14)
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
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_',ModelCollection{ivar},'_COMP-ALL.eps']);
    
end

axis off
l1=legend(ph,'-x1--','-x2--','-x3--','-Ref--');
set(l1,'Location','SouthEast','Orientation','horizontal','FontSize',14)
l1.Position = [l1.Position(1)+0.1 l1.Position(2)+0.1 l1.Position(3:4)];
xlim([Results(1).t(1) Results(1).t(end)])
ylabel('xi','FontSize',14)
xlabel('Time','FontSize',14)
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 50])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_Model_legend_COMP-ALL.eps']);

%% Show results of unforced system

ylimits = [0 11; 0 9; 0 1500];
clear ph hh
figure, hold on, box on
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
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_UNFORCED','_COMP-ALL.eps']);

%% Plot legend
% Version 1
set(gcf,'Position',[100 100 300 200])
legend off
axis on
xlim([Results(1).t(end) Results(1).t(end)+1])
ylim([10000 20000])
l1=legend(ph,'-x1--','-x2--','-x3--','-Ref--');
set(l1,'Location','SouthEast','Orientation','vertical','FontSize',14)
axis off
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 120 150])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_MPC_Comparison_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'_Model_legend_COMP-ALL_short.eps']);

% Version 2
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


%% Save results
save(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_N_',num2str(N),'_COMP-ALL.mat']),'Results')
save(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_All','_N_',num2str(N),'_COMP-ALL.mat']))