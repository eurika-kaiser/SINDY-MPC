% LOTKA-VOLTERRA system

clear all, close all, clc
figpath = '../FIGURES/';
addpath('../utils');

%% Generate Data
% Parameters: SINDy
polyorder = 3;
usesine = 0;
lambda = 0.001;      % lambda is our sparsification knob.

% Parameters: Model
a = .5;
b = .025;
d = .005;
g = .5;
n = 2;
x0=[100; 50];
dt = .01;

forcing = @(x,t) [(2*(sin(1*t)+sin(.1*t))).^2];

% Integrate
tspan=[dt:dt:100];
% u = sin(1*tspan)+sin(.1*tspan) + 2;
u = forcing(0,tspan);
N = length(tspan);
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);

%%
figure;box on
plot(t,x,'LineWidth',1)
xlabel('Time')
ylabel('Population size')
legend('Prey','Predator')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', [figpath,'EX_LOTKA_Dynamics_Noise.eps']);

%% SINDYc // CLEAN
xclean = x;
dxclean = zeros(size(xclean));
% Compute clean derivative  (just for comparison!)
for i=1:length(x)
    dxclean(i,:) = lotkacontrol(0,xclean(i,:),u(i),a,b,d,g);
end
xclean = [xclean u'];
dxclean(:,3) = 0*dxclean(:,2);
n = 3;

% SINDY+CONTROL
clear Theta Xi
Theta = poolData(xclean,n,polyorder,usesine);
m = size(Theta,2);
Xiclean = sparsifyDynamics(Theta,dxclean,lambda,n);
poolDataLIST({'x','y','u'},Xiclean,n,polyorder,usesine);

x2_std = std(xclean(:,2));
eps_vec = [0.01 0.05 0.1 0.2 0.3 0.4 0.5].*x2_std;

for i = 1:length(eps_vec)
    close all
    %% SINDYc // Noisy
    % add noise
    clear dxt xt
    eps = eps_vec(i);
    x = x + eps*randn(size(x));
    xmax = max(x);
    %  Total Variation Regularized Differentiation
    dxt(:,1) = TVRegDiff( x(:,1), 20, .00002, [], 'small', 1e12, dt, 1, 1 ); %.00002
    hold on
    plot(dxclean(:,1),'r')
    xlim([5000 7500])
    figure
    dxt(:,2) = TVRegDiff( x(:,2), 20, .00002, [], 'small', 1e12, dt, 1, 1 );
    hold on
    plot(dxclean(:,2),'r')
    xlim([5000 7500])
    
    xt(:,1) = cumsum(dxt(:,1))*dt; %xt(:,1) = xt(:,1) + x0(1);
    xt(:,2) = cumsum(dxt(:,2))*dt; %xt(:,2) = xt(:,2) + x0(2);
    xt(:,1) = xt(:,1) - (mean(xt(1000:end-1000,1)) - mean(x(1000:end-1000,1)));
    xt(:,2) = xt(:,2) - (mean(xt(1000:end-1000,2)) - mean(x(1000:end-1000,2)));
    xt = xt(500:end-501,:);
    dxt = dxt(500:end-501,:);  % trim off ends (overly conservative)
    
    close all
    %% Show
    figure,
    hold on, box on
    ccolors = get(gca,'colororder');
    plot(tspan,dxclean(:,1),'-','Color',ccolors(1,:),'LineWidth',2)
    plot(tspan,dxclean(:,2),'-','Color',ccolors(2,:),'LineWidth',2)
    plot(tspan(500:end-500),dxt(1:end,1),'-','Color',ccolors(1,:)-[0 0.3 0.09],'LineWidth',1)
    plot(tspan(500:end-500),dxt(1:end,2),'-','Color',ccolors(2,:)-[0 0.3 0.09],'LineWidth',1)
    xlabel('Time')
    ylabel('Slope of population size')
    if i == 1
        legend('Prey (clean)','Predator (clean)','Prey (TVRegDiff)','Predator (TVRegDiff)','Location','SouthEast')
    end
    ylim([-100 50])
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 600 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', [figpath,['EX_LOTKA_Dynamics_Noise_TVRegDiff_dx_',num2str(i),'.eps']]);
    print('-dpdf', [figpath,['EX_LOTKA_Dynamics_Noise_TVRegDiff_dx_',num2str(i),'.pdf']]);
    
    figure,
    hold on, box on
    ccolors = get(gca,'colororder');
    plot(tspan,x(:,1),'-','Color',ccolors(1,:)+[0.15 0.2 0.2],'LineWidth',2)
    plot(tspan,x(:,2),'-','Color',ccolors(2,:)+[0.15 0.2 0.2],'LineWidth',2)
    plot(tspan,xclean(:,1),'-','Color',ccolors(1,:),'LineWidth',2)
    plot(tspan,xclean(:,2),'-','Color',ccolors(2,:),'LineWidth',2)
    plot(tspan(500:end-500),xt(1:end,1),'-','Color',ccolors(1,:)-[0 0.3 0.09],'LineWidth',1)
    plot(tspan(500:end-500),xt(1:end,2),'-','Color',ccolors(2,:)-[0 0.3 0.09],'LineWidth',1)
    xlabel('Time')
    ylabel('Population size')
    if i == 1
        legend('Prey (corrupted)','Predator (corrupted)','Prey (clean)','Predator (clean)','Prey (TVRegDiff)','Predator (TVRegDiff)')
    end
    ylim([0 150])
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 500 200])
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', [figpath,['EX_LOTKA_Dynamics_Noise_TVRegDiff_',num2str(i),'.eps']]);
    print('-dpdf', [figpath,['EX_LOTKA_Dynamics_Noise_TVRegDiff_',num2str(i),'.pdf']]);
    %%
    xt = [xt(1:end,:) u(500:end-500)'];
    dxt = dxt(1:end,:);
    dxt(:,3) = 0*dxt(1:end,2);
    
    %% SINDY+CONTROL
    clear Theta Xi
    Theta = poolData(xt,n,polyorder,usesine);
    m = size(Theta,2);
    Xinoise = sparsifyDynamics(Theta,dxt,lambda,n);
    poolDataLIST({'x','y','u'},Xinoise,n,polyorder,usesine);
    
    %% FIGURE 1:  LORENZ for T\in[0,20]
    x0 = x(end,1:2);
    tspan_test = [100:dt:200];
    [tA,xA]=ode45(@(t,x)lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan_test,x0,options);   % true model
    [tB,xB]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xiclean(:,1:2),polyorder,usesine),tspan_test,x0,options);  % approximate // clean data
    [tC,xC]=ode45(@(t,x)sparseGalerkinControl(t,x,forcing(x,t),Xinoise(:,1:2),polyorder,usesine),tspan_test,x0,options);  % approximate // noisy data
    
    %%
    h = figure;
    subplot(2,1,1), box on
    plot(t,x(:,1),'Color',[.4 .4 .4],'LineWidth',1.5), hold on
    plot(tA,xA(:,1),'k','LineWidth',1.5), hold on
    plot(tB,xB(:,1),'r--','LineWidth',1.5)
    plot(tB,xC(:,1),'b-.','LineWidth',1.5)
    grid on
    ylim([0 200])
    % xlabel('Time','FontSize',13)
    ylabel('Prey, x_1','FontSize',13)
    set(gca,'FontSize',13)
    subplot(2,1,2), box on
    plot(t,x(:,2),'Color',[.4 .4 .4],'LineWidth',1.5), hold on
    plot(tA,xA(:,2),'k','LineWidth',1.5), hold on
    plot(tB,xB(:,2),'r--','LineWidth',1.5)
    plot(tB,xC(:,2),'b-.','LineWidth',1.5)
    if i == 1
        l1=legend('Training','Validation','SINDYc (clean)', 'SINDYc (TVRegDiff)');
        set(l1,'Location','NorthWest')
    end
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
    print(h,'-dpdf', [figpath,'EX_LOTKA_ControlValidation_Noise_',num2str(i),'.pdf'],'-r0');
    print(h,'-depsc2', [figpath,'EX_LOTKA_ControlValidation_Noise_',num2str(i),'.eps'],'-r0');
    
end