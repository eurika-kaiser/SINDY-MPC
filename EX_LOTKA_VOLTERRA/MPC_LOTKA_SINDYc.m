% MPC applied to LOTKA-VOLTERRA system using
% a SINDYc model for different prediction horizon lengths.


clear all, close all, clc
figpath = '../FIGURES/';
datapath = '../DATA/';
addpath('../utils');

%% Generate Data
% Parameters: SINDy
polyorder = 3;
usesine = 0;

% True Parameters of Lotka-Volterra model
a = .5;
b = .025;
d = .005;
g = .5;
n = 2;
x0=[60; 50];
dt = .01;

% Choose forcing function to excite system for model identification
% forcing = @(x,t) [(0.33*(sin(1*t)+sin(.1*t)))];
forcing = @(x,t) [(2*(sin(1*t)+sin(.1*t))).^2];

% Integrate excited system
tspan=[0:dt:100];
u = forcing(0,tspan);
N = length(tspan);
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);
plot(t,x,'LineWidth',1.5)
xlabel('Time')
ylabel('Population size')
legend('Prey','Predator')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', [figpath,'EX_LOTKA_Dynamics.eps']);

%% SINDYc Model Identification
% Compute Derivative
xold = x;
clear dx
x = xold;
eps = 0.1;
for i=1:length(x)
    dx(i,:) = lotkacontrol(0,x(i,:),u(i),a,b,d,g);
end
x = [xold u'];
dx(:,3) = 0*dx(:,2);
dx = dx + eps*randn(size(dx));
n = 3;

% Sparse regression
clear Theta Xi
Theta = poolData(x,n,polyorder,usesine);
m = size(Theta,2);

lambda = 0.001;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n);
poolDataLIST({'x','y','u'},Xi,n,polyorder,usesine);


%% FIGURE 1:  Lotka-Volterra // Validation
x0 = x(end,1:2);
tspan = [100 200];
[tA,xA]=ode45(@(t,x)lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);   % true model
[tB,xB]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xi(:,1:2),polyorder,usesine),tspan,x0,options);  % approximate

h = figure;
subplot(2,1,1), box on
plot(t,x(:,1),'Color',[.4 .4 .4],'LineWidth',1.5), hold on
plot(tA,xA(:,1),'k','LineWidth',1.5), hold on
plot(tB,xB(:,1),'r--','LineWidth',1.5)
grid on
ylim([0 110])
% xlabel('Time','FontSize',13)
ylabel('Prey, x_1','FontSize',13)
set(gca,'FontSize',13, 'LineWidth',1)
subplot(2,1,2), box on
plot(t,x(:,2),'Color',[.4 .4 .4],'LineWidth',1.5), hold on
plot(tA,xA(:,2),'k','LineWidth',1.5), hold on
plot(tB,xB(:,2),'r--','LineWidth',1.5)
l1=legend('Training','Validation','SINDYc');
set(l1,'Location','NorthWest')
grid on
ylim([0 60])
ylabel('Predator, x_2','FontSize',13)
set(gca,'FontSize',13, 'LineWidth',1)
xlabel('Time','FontSize',13)
set(gca,'FontSize',13)

set(h,'Units','Inches');
set(gcf,'Position',[1 1 6. 5.5])
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
print(h,'-dpdf', [figpath,'EX_LOTKA_ControlValidation.pdf'],'-r0');
print(h,'-depsc2', [figpath,'EX_LOTKA_ControlValidation.eps'],'-r0');


clear ph
figure;hold on, box on,
ccolors = get(gca,'colororder');
plot([100,100],[0 260],':','Color',[0.4,0.4,0.4],'LineWidth',1)
text(5,120,'Training', 'FontSize',12)
text(105,120,'Prediction', 'FontSize',12)
plot(t,x(:,1),'Color',[.4 .4 .4],'LineWidth',1.5); hold on
plot(t,x(:,2),'Color',[.4 .4 .4],'LineWidth',1.5)
ph(1) = plot(tA,xA(:,1),'k','LineWidth',1.5);
plot(tA,xA(:,2),'k','LineWidth',1.5)
ph(2) = plot(tB,xB(:,1),'--','Color',ccolors(1,:),'LineWidth',1.5);
ph(3) = plot(tB,xB(:,2),'--','Color',ccolors(2,:),'LineWidth',1.5);
xlabel('Time')
ylabel('Population size')
l1=legend(ph,'Validation','SINDYc','SINDYc');
axis tight
ylim([-15 260])
set(gca,'xtick',[50,100,150,200])
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-dpdf', [figpath,'EX_LOTKA_ControlValidation2']);
print('-depsc2', [figpath,'EX_LOTKA_ControlValidation2']);


%% Apply Model predictive controlto system using SINDYc model
pest.ahat = Xi(:,1:2);
pest.polyorder = polyorder;
pest.usesine = usesine;
p.a = a; % True model parameters
p.b = b;
p.d = d;
p.g = g;

% Choose prediction horizon over which the optimization is performed
% Nvec = [1,3,5,7,10,15,20,25,30,35,40,45,50];
Nvec = [2,4,6,8,9,11,12,13];

for i = 1:length(Nvec)
    Ts          = 0.1;              % Sampling time
    N           = Nvec(i);          % Control / prediction horizon (number of iterations)
    Duration    = 100;              % Run control for 100 time units
    Nvar        = 2;
    Q           = [1 0];            % State weights
    R           = 0.5;%0.5;         % Control variation du weights
    Ru = 0.5;%0.01;                 % Control weights
    B = [0; 1];                     % Control vector (which state is controlled)
    C = eye(Nvar);                  % Measurement matrix
    D = 0;                          % Feedforward (none)
    x0n=x0';%[100; 50];             % Initial condition
    uopt0 = 0;                      % Set initial control input to zero
    
    % Constraints on control optimization
    LB = [];%-100*ones(N,1);        % Lower bound of control input
    UB = [];%100*ones(N,1);         % Upper bound of control input
   
    % Reference state, which shall be achieved
    xref1 = [g/d;a/b]; % critical point
    % xref1 = [50;0];      % Reference values
    % xref2 = [50;0];
    % xref_vec = [xref1(1)*ones(size(0:Ts:10)),xref2(1)*ones(size(10+Ts:Ts:Duration));
    %             xref1(2)*ones(size(0:Ts:10)),xref2(2)*ones(size(10+Ts:Ts:Duration))];
    
    % Options for optimization routine
    options = optimoptions('fmincon','Algorithm','sqp','Display','none');
    
    % Start simulation
    fprintf('Simulation started.  It might take a while...\n')
    x        = x0n;
    Ton      = 30;      % Time when control starts
    uopt     = uopt0.*ones(N,1);
    xHistory = x;       % Stores state history
    uHistory = uopt(1); % Stores control history
    tHistory = 0;       % Stores time history
    rHistory = xref1;   % Stores reference (could be trajectory and vary with time)
    tic
    for ct = 1:(Duration/Ts)   % For each iteration: take measurements & optimize control input & apply control input
        if ct*Ts>30            % Turn control on
            if ct*Ts==Ton+Ts
                disp('Start control.')
            end
            
            % Set references
            xref = xref1;
            
            % NMPC with full-state feedback
            COSTFUN = @(u) lotkaObjectiveFCN(u,x,Ts,N,xref,uopt(1),pest,diag(Q),R,Ru);
            CONSFUN = @(u) lotkaConstraintFCN(u,x,Ts,N,pest);
            uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
            %  uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options);
            %  %use this without constraint functions CONSFUN
            
        else                    % If control is off
            uopt = uopt0.*ones(N,1);
            xref = [nan; nan];
        end
        
        % Integrate system: Apply control & Step one timestep forward
        x = rk4u(@lotkacontrol_discrete,x,uopt(1),Ts/1,1,[],p); %10, 2
        xHistory = [xHistory x];
        uHistory = [uHistory uopt(1)];
        tHistory = [tHistory tHistory(end)+Ts/1];
        rHistory = [rHistory xref];
        
    end
    fprintf('Simulation finished!\n')
    toc
    
    %% Show results
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
