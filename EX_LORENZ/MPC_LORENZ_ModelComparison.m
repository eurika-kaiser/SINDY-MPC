% LORENZ system
%118.2675
%98.2860
%2.1048e+03

clear all, close all, clc

SystemModel = 'LORENZ';

figpath = ['../../FIGURES/',SystemModel,'/'];mkdir(figpath)
datapath = ['../../DATA/',SystemModel,'/'];mkdir(datapath)
% datapath = ['/Users/ekaiser/Documents/Academia/Papers/KaKuBr_SINDYc-MPC/DATA/'];
addpath('../utils');

%% Load Models
Nvar = 3;
InputSignalTypeModel = 'sphs'; % prbs; chirp; noise; sine2; sphs; mixed

% DelayDMDc
ModelTypeDMDc = 'DMDc';
load(fullfile(datapath,['EX_',SystemModel,'_SI_',ModelTypeDMDc,'_',InputSignalTypeModel,'.mat'])) % Model: DelayDMDc , DMDc
Models.DelayDMDc = Model;

% SINDYc
load(fullfile(datapath,['EX_',SystemModel,'_SI_SINDYc','_',InputSignalTypeModel,'.mat'])) % Model
Models.SINDYc = Model;

% NARX
NARX_SUBSTRACT_MEAN = 0;
load(fullfile(datapath,['EX_',SystemModel,'_SI_NARX','_',InputSignalTypeModel,'.mat'])) % Model
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
forcing = @(x,t) (5*sin(30*t)).^3;
% tspanV =[tspan(end):dt:20];

% Reference
[tA,xA] = ode45(@(t,x)LorenzSys(t,x,forcing(x,t),p),tspanV,x(end,1:3),options);   % true model
u_valid = forcing(0,tA)';

% DelayDMDc / DMDc
% Model
iDMDc = 1; %%%% !!!!!!!!!!!!!!!!!!!! May use other one?
if Ndelay == 1
    x0      = [x(end,1:Nvar)];
    Hunew   = [u(end),u_valid(1:end-1)];
    [xB,tB] = lsim(Models.DelayDMDc.sys{iDMDc},Hunew,tspanV,[x0-[Models.DelayDMDc.xmean{iDMDc}]]');
elseif Ndelay > 1
    x0      = [x(end-Ndelay+1,1:Nvar),x(end,1:Nvar)];
    Hunew   = [ u(end-Ndelay+1:end),u_valid(1:end-Ndelay);
        u(end),u_valid(1:end-1)];
    [xB,tB] = lsim(Models.DelayDMDc.sys{1},Hunew,tspanV,[x0-[Models.DelayDMDc.xmean{1}]]');
    xB = xB(:,Nvar+1:2*Nvar); xB = xB + repmat(xmean,[size(xB,1) 1]);
end
xB = xB + repmat(Models.DelayDMDc.xmean{iDMDc},[length(tB) 1]);

% SINDYc
[tC,xC]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Models.SINDYc.Xi(:,1:Nvar),Models.SINDYc.polyorder,Models.SINDYc.usesine),tspanV,x(end,1:Nvar),options);  % approximate

% NARX
Udummy = [u(end),u_valid(1:end)];
Hdummy = zeros(Nvar,size(Udummy,2));
if NARX_SUBSTRACT_MEAN == 1
    NARX_xmean = Models.NARX.xmean;
else
    NARX_xmean = zeros(size(Models.NARX.xmean'));
end

Hdummy(:,1) = x(end,1:Nvar)'-NARX_xmean';
[Us,Ui,Si] = preparets(Models.NARX.net,con2seq(Udummy),{},con2seq(Hdummy));
xD = Models.NARX.net(Us,Ui,Si); % Predict on validation data
xD = cell2mat(xD)'; xD = [x0;xD(1:end-1,:)]; 
xD(2:end,:) = xD(2:end,:) + repmat(NARX_xmean,[size(xD,1)-1 1]);

%
clear ph
h = figure;
subplot(3,1,1), box on, hold on
plot([tB(1) tB(1)],[-100 100],'--k')
plot(t,x(:,1),'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,1),'k','LineWidth',1.5);
plot(tB,xB(:,1),'r-','LineWidth',1.5);
plot(tC,xC(:,1),'g--','LineWidth',1.5);
plot(tA(1:end),xD(:,1),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
grid on
ylim([-20 20])
% xlabel('Time','FontSize',13)
ylabel('x_1','FontSize',13)
set(gca,'FontSize',13)

subplot(3,1,2), box on, hold on
plot([tB(1) tB(1)],[-100 100],'--k')
plot(t,x(:,2),'Color',[.4 .4 .4],'LineWidth',1.5);
plot(tA,xA(:,2),'k','LineWidth',1.5);
plot(tB,xB(:,2),'r-','LineWidth',1.5);
plot(tC,xC(:,2),'g--','LineWidth',1.5);
plot(tA(1:end),xD(:,2),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
ylim([-20 20])
set(gca,'FontSize',13)
ylabel('x_2','FontSize',13)

subplot(3,1,3), box on, hold on
plot([tB(1) tB(1)],[-100 100],'--k')
ph(1) = plot(t,x(:,3),'Color',[.4 .4 .4],'LineWidth',1.5);
ph(2) = plot(tA,xA(:,3),'k','LineWidth',1.5);
ph(3) = plot(tB,xB(:,3),'r-','LineWidth',1.5);
ph(4) = plot(tC,xC(:,3),'g--','LineWidth',1.5);
ph(5) = plot(tA(1:end),xD(:,3),'-.','Color',[0.7,0.7,1],'LineWidth',1.5);
l1=legend(ph,'Training','Validation',ModelTypeDMDc,'SINDYc','NARX');
set(l1,'Location','NorthWest')
grid on
ylim([0 60])
ylabel('x_3','FontSize',13)
set(gca,'FontSize',13)
xlabel('Time','FontSize',13)


set(h,'Units','Inches');
set(gcf,'Position',[1 1 6. 5.5])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
print(h,'-painters','-depsc2',  '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'.eps'],'-r0');


%% Actuation signal
clear ph
figure,box on, hold on,
ccolors = get(gca,'colororder');
plot([tA(1),tA(1)],[-15 260],':','Color',[0.4,0.4,0.4],'LineWidth',1.5)
plot(t,u,'-k','LineWidth',1);
plot(t_valid,u_valid,'-k','LineWidth',1);
grid off
ylim([min([u u_valid])+0.05*min([u u_valid]) max([u u_valid])+0.05*max([u u_valid])])
xlim([0 t_valid(end)])
xlabel('Time')
ylabel('Input')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_Actuation.eps']);


%% MPC
% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
N  = 10;                        % Control / prediction horizon (number of iterations)
Duration = 5;                 % Run for 'Duration' time units
Ton = 0;                        % Time units when control is turned on
% Nvar = 3; % already defined
getMPCparams
          
x0n=xA(end,:)';                 % Initial condition
% xref1 = [sqrt(p.BETA*(p.RHO-1));sqrt(p.BETA*(p.RHO-1));p.RHO-1];     % critical point
xref1 = [-sqrt(p.BETA*(p.RHO-1));-sqrt(p.BETA*(p.RHO-1));p.RHO-1];
Tcontrol = tA(end);             % Time offset to combine training, prediction and control phase

% Parameters Models
ModelCollection = {ModelTypeDMDc, 'SINDYc', 'NARX'}; % 51.8924, 46.0083, 1.4662e+03
Nmodels = length(ModelCollection);

% Parameters True Model
p.SIGMA = 10;             % True system parameters
p.RHO = 28;
p.BETA = 8/3;

clear Results
Results(1:Nmodels) = struct('x',[], 'u', [], 't', [], 'xref', [], 'J', []);

Ts = 0.01;%Models.SINDYc.dt;
%%
for iM = 1:Nmodels
    select_model = ModelCollection{iM}; % DelayDMDc; DMDc; NARX ; SINDYc
    
    % Prepare variables
    Nt = (Duration/Ts)+1;
    uopt0    = 0;
    xhat     = x0n;
    uopt     = uopt0.*ones(N,1);
    xHistory = zeros(Nvar,Nt); xHistory(:,1) = xhat;
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
            %Ts = Models.DelayDMDc.dt; % Sampling time
        case 'DMDc'
            Ndmd = length(Models.DelayDMDc.sys);
            DMDc_mean = zeros(Ndmd,Nvar);
            for iDMD = 1:Ndmd
                DMDc_mean(iDMD,:) = Models.DelayDMDc.xmean{iDMD};
            end
%             [~,idx] = min( sum( abs( DMDc_mean-repmat(x0n',[Ndmd 1]) ).^2 ,2) );
            [~,idx] = min( sum( abs( DMDc_mean-repmat(xref1',[Ndmd 1]) ).^2 ,2) );
            pest.dt = Models.DelayDMDc.dt;
            pest.sys = Models.DelayDMDc.sys{idx};
            pest.xmean = DMDc_mean(idx,:)';
            pest.Nxlim = size(x,1);
            %Ts = Models.DelayDMDc.dt; % Sampling time    
        case 'SINDYc'
            pest.ahat = Models.SINDYc.Xi(:,1:Nvar);
%             pest.ahat([2,3,5],1) = [-10,10,1];
%             pest.ahat([2,3,8],2) = [28,-1,-1];
%             pest.ahat([4,7],3) = [-8/3,1];
            pest.polyorder = Models.SINDYc.polyorder;
            pest.usesine = Models.SINDYc.usesine;
            pest.dt = Models.SINDYc.dt;%*10
            %Ts = Models.SINDYc.dt; % Sampling time
            
        case 'NARX'
            pest.net = Models.NARX.net;
            pest.Ndelay = 1;%Models.NARX.Ndelay;
            pest.udelay = zeros(1,Ndelay);
            pest.xdelay = zeros(2,length(pest.udelay));
            pest.xdelay(:,1:Ndelay) = [x(end-Ndelay+1:end,1:2)]';
            pest.Nxlim = size(x,1);
            %pest.xmean = xmean';
%             pest.umean = Models.NARX.umean;
            %Ts = Models.NARX.dt; % Sampling time
    end
    
    
    % Start simulation
    fprintf('Simulation started.  It might take a while...\n')
    tic
    for ct = 1:(Duration/Ts)
        % Set references
        xref = xref1;
        
        % NMPC with full-state feedback
        COSTFUN = @(u) ObjectiveFCN_models(u,xhat,N,xref,uopt(1),pest,diag(Q),R,Ru,select_model);
        uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options);
        
        
        % Integrate system
        xhat = rk4u(@lorenzcontrol_discrete,xhat,uopt(1),Ts/10,10,[],p); %10, 2
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
%         elseif strcmp(select_model, 'DMDc')
%             [~,idx] = min( sum( abs( DMDc_mean-repmat(xhat',[Ndmd 1]) ).^2 ,2) );
%             pest.sys = Models.DelayDMDc.sys{idx};
%             pest.xmean = DMDc_mean(idx,:)';
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

%N=10: 120.9128, 96.8967, 2.1755e+03

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

%% Unforced
[tU,xU] = ode45(@(t,x)LorenzSys(t,x,0,p),tHistory,xA(end,1:3),options);   % true model

%% Show cost
clear ph
figure, hold on, box on
ph(1) = plot(Results(1).t,cumsum(Results(1).J),'-r','LineWidth',2);
ph(3) = plot(Results(1).t,cumsum(Results(3).J),'-.','Color',[0.7,0.7,1],'LineWidth',2);
ph(2) = plot(Results(1).t,cumsum(Results(2).J),'--g','LineWidth',2);

l1=legend(ph,ModelCollection);
set(l1,'Location','SouthEast')
xlim([Results(1).t(1) Results(1).t(end)])
% ylim([0 60])
ylabel('Cost','FontSize',14)
xlabel('Time','FontSize',14)
% set(gca,'yscale','log')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

print('-depsc2', '-loose','-cmyk', [figpath,'EX_',SystemModel,'_SI_Comparison_Validation_',InputSignalTypeModel,'_',InputSignalType,'_Cost','_N_',num2str(N),'.eps']);
    

%% Save results
save(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_N_',num2str(N),'.mat']),'Results')
save(fullfile(datapath,['MPC_',SystemModel,'_Results_AllModels_',InputSignalTypeModel,'_',InputSignalType,'_All','_N_',num2str(N),'.mat']))