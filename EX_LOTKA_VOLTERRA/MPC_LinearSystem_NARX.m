% LOTKA-VOLTERRA system

clear all, close all, clc
figpath = '../../FIGURES/';
datapath = '../../DATA/';
addpath('../utils');

%% Generate Data
sys = 'LinSys';

% System
dt = .01;
omega = 0.1;
sigma = 0.04;
LinearSystem = @(t,x,u) [-sigma omega; -omega -sigma]*x + [1;u];
Pf = 10; K = 10; A = 5;
forcing = @(x,t) A*sphs(Pf,K,t); %(t<50)*A*sphs(Pf,K,t);

% Integrate
x0 = [1;1];
tspan=[dt:dt:100];
N = length(tspan);
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,2));
[t,x]=ode45(@(t,x) LinearSystem(t,x,forcing(x,t)),tspan,x0,options);


plot(x(:,1),x(:,2),'LineWidth',1.5)
xlabel('x')
ylabel('y')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', [figpath,'EX_',sys,'_Dynamics.eps']);

u = forcing(0,tspan);
figure,plot(u)

x = x';
Nt = length(tspan);
Ndelay = 1;
%% NARX: Training
rng(2,'twister')

% prepare training data
yt = con2seq(x);
yi = con2seq(u); 

% prepare validation data
%yv = con2seq(xA');

% Neural network
stateDelays = 1;
inputDelays = 1;
hiddenSizes = [10];

% Nonlinear autoregressive neural network
net = narxnet(inputDelays,stateDelays, hiddenSizes);

% Training parameters %nnstart
net.trainFcn = 'trainlm';%trainbr; trainlm; trainscg
net.trainParam.min_grad = 1e-10;
net.trainParam.showCommandLine = 1;
% net.trainParam.epochs = 1000;
% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 15/100;
% net.performFcn = 'mse';  % Mean squared error

% Prepares training data (shifting, copying feedback targets into inputs as needed, etc.)
[Us,Ui,Si,Ss] = preparets(net,yi,{},yt); %yt

% Train net with prepared training data in open-loop
tic
net = train(net,Us,Ss,Ui,Si);
toc

% Plots
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotregression(t,y)
%figure, plotresponse(t,y)
%figure, ploterrcorr(e)
%figure, plotinerrcorr(x,e)

% Close loop for recursive prediction
netc = closeloop(net);

%% NARX: Prediction
% Prepare validation data / Get initial state from training data
[Us,Ui,Si,So] = preparets(netc,yi,{},yt); 

% Predict on validation data
predict = netc(Us,Ui,Si);
xNARX = cell2mat(predict)';

% Error
e = cell2mat(gsubtract(So,predict)); 

%
figure;
plot(tspan(max(stateDelays):Nt-1),xNARX','-','LineWidth',1,'Color','k');%0.7*ones(1,3))
grid on, hold on

%% Show 
clear ph
figure,box on,
ccolors = get(gca,'colororder');
ph(1) = plot(tspan,x(1,:),'-','Color',ccolors(1,:),'LineWidth',1); hold on
ph(2) = plot(tspan,x(2,:),'-','Color',ccolors(2,:),'LineWidth',1);
ph(3) = plot(tspan(Ndelay+1:Nt),xNARX(:,1),'--','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',2);
ph(4) = plot(tspan(Ndelay+1:Nt),xNARX(:,2),'--','Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',2);
xlim([0 100])
xlabel('Time')
ylabel('x_i')
legend(ph([1,3]),'True','NARX')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', [figpath,'EX_',sys,'_SI_NARX.eps']);

%% MPC
% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
Ts = 0.1;                   % Sampling time
N  = 10;                     % Control / prediction horizon (number of iterations)
Duration = 10;             % Run for 'Duration' time units
Nt = (Duration/Ts)+1;       %
Ton = 0;                    % Time units when control is turned on
Nvar = 2;
Q = [1 0];                  % State weights
R = 0.5;%0.5;                    % du weights
Ru = 0.5;%0.5;                   % u weights
B = [0; 1];
C = eye(Nvar);
D = 0;
LB = [];%-100*ones(N,1);    % Lower bound of control input
UB = [];%100*ones(N,1);     % Upper bound of control input
x0n=x0(3:4)';                     % Initial condition
xref1 = [g/d;a/b];          % critical point

% Parameters Models
ModelCollection = {'DelayDMDc', 'SINDYc', 'NARX'};

% Parameters True Model
p.a = a;
p.b = b;
p.d = d;
p.g = g;


%% MPC
pest.ahat = Xi(:,1:2);
pest.polyorder = polyorder;
pest.usesine = usesine;
p.a = a;
p.b = b;
p.d = d;
p.g = g;

%%
% Nvec = [1,3,5,7,10,15,20,25,30,35,40,45,50];
Nvec = [2,4,6,8,9,11,12,13];

for i = 1:length(Nvec);
    Ts          = 0.1;           % Sampling time
    N           = Nvec(i);             % Control / prediction horizon (number of iterations)
    Duration    = 100;     % Run for 100 time units
    Nvar        = 2;
    Q           = [1 0];          % State weights
    R           = 0.5;%0.5;            % du weights
    Ru = 0.5;%0.01;     % Control weights
    B = [0; 1];
    C = eye(Nvar);
    D = 0;
    LB = [];%-100*ones(N,1);    % Lower bound of control input
    UB = [];%100*ones(N,1);     % Upper bound of control input
    x0n=x0';%[100; 50];       % Initial condition
    uopt0 = 0;
    
    
    xref1 = [g/d;a/b]; % critical point
    % xref1 = [50;0];      % Reference values
    % xref2 = [50;0];
    % xref_vec = [xref1(1)*ones(size(0:Ts:10)),xref2(1)*ones(size(10+Ts:Ts:Duration));
    %             xref1(2)*ones(size(0:Ts:10)),xref2(2)*ones(size(10+Ts:Ts:Duration))];
    
    options = optimoptions('fmincon','Algorithm','sqp','Display','none');
    
    % Start simulation
    fprintf('Simulation started.  It might take a while...\n')
    x = x0n;
    Ton = 30;
    uopt = uopt0.*ones(N,1);
    xHistory = x;
    uHistory = uopt(1);
    tHistory = 0;
    rHistory = xref1;
    tic
    for ct = 1:(Duration/Ts)
        if ct*Ts>30 % turn control on
            if ct*Ts==Ton+Ts
                disp('Start control.')
            end
            % Set references.
            xref = xref1;
            
            % NMPC with full-state feedback
            COSTFUN = @(u) lotkaObjectiveFCN(u,x,Ts,N,xref,uopt(1),pest,diag(Q),R,Ru);
            CONSFUN = @(u) lotkaConstraintFCN(u,x,Ts,N,pest);
            uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
            %         uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options);
            %     keyboard
        else
            uopt = uopt0.*ones(N,1);
            xref = [nan; nan];
        end
        
        % Integrate system
        x = rk4u(@lotkacontrol_discrete,x,uopt(1),Ts/1,1,[],p); %10, 2
        xHistory = [xHistory x];
        uHistory = [uHistory uopt(1)];
        tHistory = [tHistory tHistory(end)+Ts/1];
        rHistory = [rHistory xref];
        
    end
    fprintf('Simulation finished!\n')
    toc
    %%
    clear ph
    
    figure;hold on, box on,
    ccolors = get(gca,'colororder');
    plot([Ton+tspan(1),Ton+tspan(1)],[-15 260],':','Color',[0.4,0.4,0.4],'LineWidth',1)
    text(31+tspan(1),210,'Control', 'FontSize',12)
    text(31+tspan(1),190,'turned on', 'FontSize',12)
    plot(tHistory+tspan(1),zeros(length(tHistory),1),'-k','LineWidth',0.5)
    plot(tHistory+tspan(1),xref1(1)*ones(length(tHistory),1),'--','Color',ccolors(1,:),'LineWidth',1)
    plot(tHistory+tspan(1),xref1(2)*ones(length(tHistory),1),'--','Color',ccolors(2,:),'LineWidth',1)
    ph(1) = plot(tHistory+tspan(1),xHistory(1,:),'-','Color',ccolors(1,:),'LineWidth',1.5);
    ph(2) = plot(tHistory+tspan(1),xHistory(2,:),'-','Color',ccolors(2,:),'LineWidth',1.5);
    ph(3) = plot(tHistory+tspan(1),uHistory,'-k','LineWidth',1.5);
    xlabel('Time')
    ylabel('Population size')
    legend(ph,'Prey','Predator','Control')
    axis tight
    ylim([-15 260])
    xlim([100,200.0001])
    %ylim([min(uHistory)-5 260])
    set(gca,'xtick',[50,100,150,200])
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', [figpath,'EX_LOTKA_DynamicsControlled_cnstrnd_N',num2str(N),'.eps']);
    
    %% Save Results
    Results.t = tHistory;
    Results.x = xHistory;
    Results.u = uHistory;
    Results.J = evalObjectiveFCN(uHistory,xHistory,rHistory,diag(Q),R,Ru);
    
    save(fullfile(datapath,['EX_LOTKA_MPC_SINDYc_N',num2str(N),'.mat']),'Results')
    
end
return
%% Timestep
Ts = 0.1;
Lotka = @(t,x) [a*x(1)-b*x(1)*x(2);
    -g*x(2)+d*x(1)*x(2)];
Duration = 100;
xnl = x0;
for i = 1:Duration/Ts
    xnlnew = rk4u(@lotkacontrol_discrete,xnl(:,end),0,Ts,1,[],p);
    xnl = [xnl, xnlnew];
end

[t,xt] = ode45(Lotka,[1:Duration/Ts].*Ts,x0);

figure,
plot([0:Duration/Ts].*Ts,xnl','-k'), hold on, plot([0:Duration/Ts-1].*Ts,xt,'-r');
%% Linearization
xref1 = [g/d,a/b];
xeq = x0';
Ts = 0.1;
Jacobian = @(x)[a-b*x(2) -b*x(1);
    d*x(2)   d*x(1)-g];
Lotka = @(x) [a*x(1)-b*x(1)*x(2);
    -g*x(2)+d*x(1)*x(2)];
dfdx = Lotka(xeq);
A = Jacobian(xeq);

sys = ss(A,B,eye(2),0);
plant = c2d(sys,Ts);

xt = x0;
xnl = x0;
for i = 1:1000
    xtnew = plant.A*(xt(:,end));
    xt = [xt, xtnew];
    xnlnew = rk4u(@lotkacontrol_discrete,xnl(:,end),0,Ts,1,[],p);
    xnl = [xnl, xnlnew];
end
figure,
plot(xnl(:,:)','-k'), hold on, plot(xt'+repmat(xeq,[size(xt',1),1]),'-r');
return