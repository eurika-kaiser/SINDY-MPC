% MPC applied to HIV system using a SINDYc model.
% Function for TESTING

clear all, close all, clc
figpath = '../FIGURES/HIV/'; mkdir(figpath)
datapath = '../DATA/HIV/'; mkdir(datapath)
addpath('../utils');

SystemModel = 'HIV';
ModelName = 'SINDYc';

%% Load Model
Nvar = 5;
InputSignalTypeModel = 'prbs'; % prbs; chirp; noise; sine2; sphs; mixed
load(fullfile(datapath,['EX_',SystemModel,'_SI_SINDYc','_',InputSignalTypeModel,'.mat'])) 

%% TRUE SYSTEM PARAMETERS
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


%% Apply Model predictive control to system using SINDYc model
select_model = 'SINDYc';
pest.ahat = Model.Xi(:,1:Nvar); % Use Xi0(:,1:Nvar) to execute true model parameters
pest.polyorder = Model.polyorder;
pest.usesine = Model.usesine;
pest.dt = 1/12; %Model.dt; % Don't need to use the same time step model was trained on


% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
Nweeks = 50;
Duration = Nweeks*7;             % Run for 'Duration' time units
Ton = 0;                         % Time units when control is turned on   
x0n = [10, 0.1, 0.1, 0.1, 0.1]'; % Initial condition
Ts  = pest.dt;                   % Sampling time [1/hr]
Tcontrol = 0;
getMPCparams    

% Reference state, which shall be achieved
xref = zeros(1,5);
xref(2) = ((c2*(lambda1-d*q)-b2*alpha1) - sqrt((c2*(lambda1-d*q)-b2*alpha1)^2 - 4*alpha1*c2*q*d*b2))/(2*alpha1*c2*q);
xref(1) = lambda1/(d+alpha1*xref(2));
xref(4) = 0;
xref(5) = (xref(2)*c2*(alpha1*q-a) + b2*alpha1)/(c2*p2*xref(2));
xref(3) = h*xref(5)/(c2*q*xref(2));


% Initialize variables
Nt = (Duration/Ts)+1;
uopt0    = 0;
xhat     = x0n;
uopt     = uopt0.*ones(Nu,1);
xHistory = zeros(Nvar,Nt); xHistory(:,1) = xhat;
uHistory = zeros(1,Nt); uHistory(1)   = uopt(1);
tHistory = zeros(1,Nt); tHistory(1)   = Tcontrol;
rHistory = zeros(Nvar,Nt);

%%
% Start simulation
fprintf('Simulation started.  It might take a while...\n')
tic
for ct = 1:(Duration/Ts)   % For each iteration: take measurements & optimize control input & apply control input
    
    if mod(ct*Ts,7) == 0 % Update once a week
        % NMPC with full-state feedback
        COSTFUN = @(u) ObjectiveFCN_models(u,xhat,N,Nu,xref,uHistory(:,ct),pest,diag(Q),R,Ru,select_model);
        CONSFUN = @(u) ConstraintFCN_models(u,uHistory(:,ct),xhat,N,LBo,UBo,LBdu,UBdu,pest,select_model);
        uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
    end
    
    % Integrate system
    xhat = rk4u(@HIVsys_ZURAKOWSKI,xhat,uopt(1),Ts/2,2,[],0); 
    xHistory(:,ct+1) = xhat;
    uHistory(:,ct+1) = uopt(1);
    tHistory(:,ct+1) = ct*Ts+Tcontrol; 
    rHistory(:,ct+1) = xref;
    
    if mod(ct,1000) == 0
        disp(['PROGRESS: ',num2str(100*ct/(Duration/Ts)),'%'])
    end
    
    
end
fprintf('Simulation finished!\n')
toc

%%
clear ph
figure;box on,hold on
ccolors = get(gca,'colororder');
plot(tHistory(1:ct)/7,xref(1)*ones(length(tHistory(1:ct)),1)./max(xHistory(1,1:ct)),'--','Color',ccolors(1,:),'LineWidth',1), hold on
plot(tHistory(1:ct)/7,xref(2)*ones(length(tHistory(1:ct)),1)./max(xHistory(2,1:ct)),'--','Color',ccolors(2,:),'LineWidth',1)
plot(tHistory(1:ct)/7,xref(3)*ones(length(tHistory(1:ct)),1)./max(xHistory(3,1:ct)),'--','Color',ccolors(3,:),'LineWidth',1)
ph(1) = plot(tHistory(1:ct)/7,xHistory(1,1:ct)./max(xHistory(1,1:ct)),'-','Color',ccolors(1,:),'LineWidth',1.5);
ph(2) = plot(tHistory(1:ct)/7,xHistory(2,1:ct)./max(xHistory(2,1:ct)),'-','Color',ccolors(2,:),'LineWidth',1.5);
ph(3) = plot(tHistory(1:ct)/7,xHistory(3,1:ct)./max(xHistory(3,1:ct)),'-','Color',ccolors(3,:),'LineWidth',1.5);
ph(4) = plot(tHistory(1:ct)/7,uHistory(1:ct),'-k','LineWidth',1.5);
ylim([0 1])
legend(ph,'x1','x2','x3','u')
%% Show results
clear ph

figure;hold on, box on,
ccolors = get(gca,'colororder');
plot(tHistory,xref1(1)*ones(length(tHistory),1),'--','Color',ccolors(1,:),'LineWidth',1)
plot(tHistory,xref1(2)*ones(length(tHistory),1),'--','Color',ccolors(2,:),'LineWidth',1)
plot(tHistory,xref1(3)*ones(length(tHistory),1),'--','Color',ccolors(3,:),'LineWidth',1)
ph(1) = plot(tHistory,xHistory(1,:),'-','Color',ccolors(1,:),'LineWidth',1.5);
ph(2) = plot(tHistory,xHistory(2,:),'-','Color',ccolors(2,:),'LineWidth',1.5);
ph(3) = plot(tHistory,xHistory(3,:),'-','Color',ccolors(3,:),'LineWidth',1.5);
ph(4) = plot(tHistory,uHistory,'-k','LineWidth',1.5);
xlabel('Time')
ylabel('xi')
legend(ph,'x1','x2','x3','Control')
axis tight
set(gca,'xtick',[50,100,150,200])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

%% Save Results
Results.t = tHistory;
Results.x = xHistory;
Results.u = uHistory;
Results.J = evalObjectiveFCN(uHistory,xHistory,rHistory,diag(Q),R,Ru);

