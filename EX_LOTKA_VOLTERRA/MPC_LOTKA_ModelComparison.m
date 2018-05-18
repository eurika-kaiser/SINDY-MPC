% LOTKA-VOLTERRA system

clear all, close all, clc
figpath = '../../FIGURES/';
datapath = '../../DATA/';
% datapath = ['/Users/ekaiser/Documents/Academia/Papers/KaKuBr_SINDYc-MPC/DATA/'];
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
NARX_SUBSTRACT_MEAN = 0;
load(fullfile(datapath,['EX_LOTKA_SI_NARX','_',InputSignalTypeModel,'.mat'])) % Model
Models.NARX = Model;

%% Get training data
ONLY_TRAINING_LENGTH = 1;
InputSignalType = InputSignalTypeModel;% prbs; chirp; noise; sine2; sphs; mixed
Ndelay = Models.DelayDMDc.Ndelay;
getTrainingData

%% Get validation data
InputSignalType = 'sine2';% prbs; chirp; noise; sine2; sphs; mixed
Ndelay = Models.DelayDMDc.Ndelay;
getValidationData


%% FIGURE 1:  // Model comparison + Validation
% Forcing
% forcing = @(x,t) [(2*(sin(1*t)+sin(.1*t))).^2]; %0.33
% forcing = @(x,t) [(1.8*(sin(1.1*t)+sin(.2*t))).^2]; %0.33
% forcing = @(x,t) Amp*chirp(t,[],max(tspan),1.).^2;

% Reference
%x0      = [x(end-Ndelay+1,1:2),x(end,1:2)];
% tspanV   = [100:dt:200];
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
% if strcmp(ModelType,'DelayDMDc')
%     unew    = forcing(0,tspan);
%     Hunew   = [ u(end-Ndelay+1:end),unew(1:end-Ndelay);
%         u(end),unew(1:end-1)];
%
%     [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,tspan,[x0-[xmean,xmean]]');
%     xB = xB(:,3:4); xB = xB + repmat(xmean,[size(xB,1) 1]);
% elseif strcmp(ModelType,'DMDc')
%     [xDMDc,~] = lsim(sysmodel_DMDc,Hu,tspan(1:Nt),Hx(:,1));
%     xDMDc = xDMDc(:,end-1:end);
%     xDMDc = xDMDc + repmat(xmean,[Nt 1]);
% end
% SINDYc
[tC,xC]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Models.SINDYc.Xi(:,1:2),Models.SINDYc.polyorder,Models.SINDYc.usesine),tspanV,x(end,1:2),options);  % approximate

% NARX
% Udummy = [u(end-Ndelay+1:end),unew(1:end)];
% Hdummy = zeros(2,size(Udummy,2));
% Hdummy(:,1:Ndelay) = x(end-Ndelay+1:end,1:2)'- repmat(xmean',[1 Ndelay]);
Udummy = [u(end),u_valid(1:end)];%-Models.NARX.umean;
Hdummy = zeros(2,size(Udummy,2));
if NARX_SUBSTRACT_MEAN == 1
    NARX_xmean = Models.NARX.xmean';
else
    NARX_xmean = zeros(size(Models.NARX.xmean'));
end

Hdummy(:,1) = x(end,1:2)'-NARX_xmean;
[Us,Ui,Si] = preparets(Models.NARX.net,con2seq(Udummy),{},con2seq(Hdummy));
xD = Models.NARX.net(Us,Ui,Si); % Predict on validation data
xD = cell2mat(xD)'; xD = [x0;xD(1:end-1,:)]; 
xD(2:end,:) = xD(2:end,:) + repmat(NARX_xmean',[size(xD,1)-1 1]);

% Error
%e = cell2mat(gsubtract(So,predict));

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
% text(tB(1)-50,max([u uv])-0.1*max([u uv]),'Training', 'FontSize',12,'Color','r')
% text(10+tB(1),max([u uv])-0.1*max([u uv]),'Validation', 'FontSize',12,'Color','r')
% text(5+tB(1),210,'turned on', 'FontSize',12)
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

% x0n=x0(3:4)';                   % Initial condition
x0n=xA(end,:)';
xref1 = [g/d;a/b];              % critical point
Tcontrol = tA(end);             % Time offset to combine training, prediction and control phase

% Parameters Models
ModelCollection = {ModelTypeDMDc, 'SINDYc', 'NARX'}; % 51.8924, 46.0083, 1.4662e+03
Nmodels = length(ModelCollection);

% Parameters True Model
p.a = a;
p.b = b;
p.d = d;
p.g = g;

clear Results
Results(1:Nmodels) = struct('x',[], 'u', [], 't', [], 'xref', [], 'J', []);

Ts = Models.SINDYc.dt;
%
for iM = 1:Nmodels
    select_model = ModelCollection{iM}; % DelayDMDc; DMDc; NARX ; SINDYc
    
    % Prepare variables
    Nt = (Duration/Ts)+1;
    uopt0    = 0;
    xhat     = x0n;
    uopt     = uopt0.*ones(N,1);
    xHistory = zeros(2,Nt); xHistory(:,1) = xhat;
    uHistory = zeros(1,Nt); uHistory(1)   = uopt(1);
    tHistory = zeros(1,Nt); tHistory(1)   = Tcontrol;
    
    % Parameters Model
    switch select_model
        case 'DelayDMDc'
            pest.dt = Models.DelayDMDc.dt;
            pest.sys = Models.DelayDMDc.sys;
            pest.xmean = xmean';
            pest.udelay = zeros(1,N);
            %         pest.xdelay = [x(end-Ndelay+1:end-Ndelay+N,1:2)]';
            pest.xdelay = [x(end-Ndelay+1,1:2)]';
            pest.Nxlim = size(x,1);
            Ts = Models.DelayDMDc.dt; % Sampling time
        case 'DMDc'
            pest.dt = Models.DelayDMDc.dt;
            pest.sys = Models.DelayDMDc.sys;
            pest.xmean = xmean';
            pest.Nxlim = size(x,1);
            Ts = Models.DelayDMDc.dt; % Sampling time    
        case 'SINDYc'
            pest.ahat = Models.SINDYc.Xi(:,1:2);
            pest.polyorder = Models.SINDYc.polyorder;
            pest.usesine = Models.SINDYc.usesine;
            pest.dt = Models.SINDYc.dt;
            Ts = Models.SINDYc.dt; % Sampling time
            
        case 'NARX'
            pest.net = Models.NARX.net;
            pest.Ndelay = 1;%Models.NARX.Ndelay;
            pest.udelay = zeros(1,Ndelay);
            pest.xdelay = zeros(2,length(pest.udelay));
            pest.xdelay(:,1:Ndelay) = [x(end-Ndelay+1:end,1:2)]';
            pest.Nxlim = size(x,1);
            pest.xmean = xmean';
%             pest.umean = Models.NARX.umean;
            Ts = Models.NARX.dt; % Sampling time
    end
    
    
    % Start simulation
    fprintf('Simulation started.  It might take a while...\n')
    tic
    for ct = 1:(Duration/Ts)
        % Set references
        xref = xref1;
        
        % NMPC with full-state feedback
        COSTFUN = @(u) lotkaObjectiveFCN_models(u,xhat,N,xref,uopt(1),pest,diag(Q),R,Ru,select_model);
        CONSFUN = @(u) lotkaConstraintFCN_models(u,xhat,N,pest,select_model);
        uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
        %   uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options);
        
        
        % Integrate system
        xhat = rk4u(@lotkacontrol_discrete,xhat,uopt(1),Ts/1,1,[],p); %10, 2
        xHistory(:,ct+1) = xhat;
        uHistory(:,ct+1) = uopt(1);
        tHistory(:,ct+1) = ct*Ts+Tcontrol;
        
        % Update states/inputs/parameters
        %     if strcmp(select_model,'DelayDMDc')
        %         if pest.Nxlim-Ndelay+N+ct<=pest.Nxlim
        %             pest.udelay = zeros(1,N);
        %             %pest.xdelay = [x(end-Ndelay+1+ct:end-Ndelay+N+ct,1:2)]'; %[x(end-Ndelay+1+ct,1:2),xhat]';
        %             pest.xdelay = [x(end-Ndelay+1+ct,1:2)]';
        %         elseif pest.Nxlim-Ndelay+1+ct<pest.Nxlim && pest.Nxlim-Ndelay+N+ct>pest.Nxlim
        %             idiff = pest.Nxlim-Ndelay+N+ct - pest.Nxlim;
        %             pest.udelay = [zeros(1,N-idiff+1),uHistory(1:ct-(Ndelay-N))];
        %             %pest.xdelay = [x(end-Ndelay+1+ct:end,1:2),xHistory(1:2,ct-(Ndelay-N))']';
        %             pest.xdelay = [x(end-Ndelay+1+ct,1:2)]';
        %         else
        %             idiff = pest.Nxlim-Ndelay+N+ct - pest.Nxlim;
        %             pest.udelay = [uHistory(ct-Ndelay+1:ct-(Ndelay-N))];
        %             %pest.xdelay = [xHistory(1:2,ct-Ndelay+1:ct-(Ndelay-N))]';
        %             pest.xdelay = [xHistory(1:2,ct-Ndelay+1)']';
        %         end
        %         if size(pest.xdelay,2)~=1
        %             disp('Error in size.')
        %             return
        %         end
        %         if length(pest.udelay)~=N
        %             disp('Error in size.')
        %             return
        %         end
        %     if strcmp(select_model,'NARX') || strcmp(select_model,'DelayDMDc')
        if strcmp(select_model,'DelayDMDc')
            if  pest.Nxlim-Ndelay+1+ct<=pest.Nxlim
                pest.udelay = [zeros(1,Ndelay-ct),uHistory(1:ct)];
                pest.xdelay = zeros(2,length(pest.udelay));
                pest.xdelay(:,1:Ndelay) = [x(end-Ndelay+1+ct:end,1:2);xHistory(1:2,1:ct)']';
            else
                pest.udelay = [uHistory(ct-Ndelay+1:ct)];
                pest.xdelay = zeros(2,length(pest.udelay));
                pest.xdelay = [xHistory(1:2,ct-Ndelay+1:ct)];
            end
            if size(pest.xdelay,2)~=Ndelay
                disp('Error in size.')
                return
            end
            if length(pest.udelay)~=Ndelay
                disp('Error in size.')
                return
            end
        end
    end
    fprintf('Simulation finished!\n')
    %
    tElapsed = toc
    
    Results(iM).eval_time = tElapsed;
    Results(iM).xref = repmat(xref,[1 length(tHistory)]);
    Results(iM).x = xHistory;
    Results(iM).u = uHistory;
    Results(iM).t = tHistory;
    Results(iM).J = evalObjectiveFCN(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);
    
    VIZ_SI_Validation_MPC
end

%N=10: 69.1819, 74.5837, 1.9040e+03
%N=5: 39.4120, 30.8995, 1.1190e+03

%% Unforced
[tU,xU] = ode45(@(t,x)lotkacontrol(t,x,0,a,b,d,g),tHistory,xA(end,1:2),options);   % true model

%% Show cost
clear ph
figure, hold on, box on
ph(1) = plot(Results(1).t,cumsum(Results(1).J),'-r','LineWidth',2);
ph(2) = plot(Results(1).t,cumsum(Results(2).J),'--g','LineWidth',2);
ph(3) = plot(Results(1).t,cumsum(Results(3).J),'-.','Color',[0.7,0.7,1],'LineWidth',2);

l1=legend(ph,ModelCollection);
set(l1,'Location','SouthEast')
% ylim([0 60])
ylabel('Cost','FontSize',14)
xlabel('Time','FontSize',14)
% set(gca,'yscale','log')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

print('-depsc2', '-loose','-cmyk', [figpath,'EX_LOTKA_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'.eps']);
    
%%
% xA = Results(2).x'; tB = tA; tC = tA; tD = tA; 
% tA = Results(2).t;
% 
% % DMDc
% x0      = [x0n'];
% Hunew   = [Results(1).u]';
% [xB,tB] = lsim(Models.DelayDMDc.sys,Hunew,Results(1).t,[x0-[xmean]]');
% xB = xB + repmat(xmean,[length(tB) 1]);
% 
% % SINDYc
% pest.ahat = Models.SINDYc.Xi(:,1:2);
% pest.polyorder = Models.SINDYc.polyorder;
% pest.usesine = Models.SINDYc.usesine;
% pest.dt = Models.SINDYc.dt;
% Ns = size(x0n,1);
% xC = zeros(Ns,length(Results(2).u)); xC(:,1) = x0n;
% for ct=1:length(Results(2).u)-1
%     % Obtain plant state at next prediction step.
%     xC(:,ct+1) = rk4u(@sparseGalerkinControl_Discrete,xC(:,ct),Results(2).u(ct),dt,1,[],pest);
% end
% xC = xC';
% %xk = xk(:,2:N+1);
%     
% % NARX
% Udummy = [0,Results(3).u];
% Hdummy = zeros(2,size(Udummy,2));
% Hdummy(:,1) = x0n;
% [Us,Ui,Si] = preparets(Models.NARX.net,con2seq(Udummy),{},con2seq(Hdummy));
% xD = Models.NARX.net(Us,Ui,Si); 
% xD = cell2mat(xD)'; xD = [x0n';xD(1:end-1,:)]; 
% 
% %
% clear ph
% h = figure;
% subplot(2,1,1), box on, hold on
% plot([tB(1) tB(1)],[0 150],'--k')
% plot(tA,xA(:,1),'k','LineWidth',1.5);
% plot(tB,xB(:,1),'r-','LineWidth',1.5);
% plot(tC,xC(:,1),'g--','LineWidth',1.5);
% plot(tA(1:end),xD(:,1),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
% grid on
% ylim([0 150])
% % xlabel('Time','FontSize',13)
% ylabel('Prey, x_1','FontSize',13)
% set(gca,'FontSize',13)
% subplot(2,1,2), box on, hold on
% plot([tB(1) tB(1)],[0 60],'--k')
% ph(1) = plot(t,x(:,2),'Color',[.4 .4 .4],'LineWidth',1.5);
% ph(2) = plot(tA,xA(:,2),'k','LineWidth',1.5);
% ph(3) = plot(tB,xB(:,2),'r-','LineWidth',1.5);
% ph(4) = plot(tC,xC(:,2),'g--','LineWidth',1.5);
% % ph(5) = plot(tB,xD(:,2),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
% ph(5) = plot(tA(1:end),xD(:,2),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
% l1=legend(ph,'Training','Validation',ModelTypeDMDc,'SINDYc','NARX');
% set(l1,'Location','NorthWest')
% grid on
% ylim([0 60]), xlim([200 220])
% ylabel('Predator, x_2','FontSize',13)
% set(gca,'FontSize',13)
% xlabel('Time','FontSize',13)
% set(gca,'FontSize',13)
% 
% set(h,'Units','Inches');
% set(gcf,'Position',[1 1 6. 5.5])
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])

%% Save results
save(fullfile(datapath,['MPC_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_N_',num2str(N),'.mat']),'Results')
save(fullfile(datapath,['MPC_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_All','_N_',num2str(N),'.mat']))