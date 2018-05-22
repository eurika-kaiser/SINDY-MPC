% LOTKA-VOLTERRA system
% Run MPC for various prediction horizon lengths

clear all, close all, clc
figpath = '../../FIGURES/';
addpath('../utils');

%% Generate Data
% Parameters: SINDy
polyorder = 3;
usesine = 0;

% Parameters: Model
a = .5;
b = .025;
d = .005;
g = .5;
n = 2;
x0=[60,50];
xref = [g/d;a/b]; % critical point

forcing = @(x,t) [(2*(sin(1*t)+sin(.1*t))).^2]; %0.33

% Integrate
dt = 0.01;
tspan=[dt:dt:100];
% u = sin(1*tspan)+sin(.1*tspan) + 2;
u = forcing(0,tspan);
N = length(tspan);
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);
plot(t,x,'LineWidth',1)
xlabel('Time')
ylabel('Population size')
legend('Prey','Predator')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
% print('-depsc2', [figpath,'EX_LOTKA_Dynamics.eps']);

xold = x;

%% DMDc: Bunknown
T = length(tspan);
% No time delay
% Gamma = u(1:T-1);
% X1 = (xold(1:T-1,:)-repmat(mean(xold(5000:T,:)),[T-1 1]))';
% X2 = (xold(2:T,:)-repmat(mean(xold(5000:T,:)),[T-1 1]))';

% % Step 1: SVD of Input space
% Omega = [X1;Gamma];
% [U1,S1,V1] = svd(Omega,'econ');
% 
% % Step 2: SVD of Output space (if projection onto subspace)
% [U2,S2,V2] = svd(X2,'econ');
% 
% % Step 3: Estimate A, B
% Atilde = U2'*X2*V1*S1^(-1)*U1(1:2,:)'*U2; 
% Btilde = U2'*X2*V1*S1^(-1)*U1(3,:)'; 
% 
% % Prediction
% sys = ss(Atilde, Btilde, eye(2),0,dt);
% [xDMDc,t] = lsim(sys,u,tspan,x0-mean(xold(1:T,:))');

% With time delay

delay_vec = [1:50];
l2err = zeros(length(delay_vec),1);
for i = 1:length(delay_vec)
    Ndelay = delay_vec(i);
    
    
    X = xold - repmat(mean(xold(5000:T,:)),[T 1]);
    Hx = getHankelMatrix_MV(X,Ndelay); Hx = Hx([1:2,2*Ndelay-1:2*Ndelay],:);
    Hu = getHankelMatrix_MV(u',Ndelay); Hu = Hu([1,Ndelay],:);
    [sysmodel_DMDc,U,Up] = DelayDMDc_MV(Hx,Hu,4,4,dt,4,2,2);

    % Prediction
    Nt = size(Hu,2);
    [xDMDc,t] = lsim(sysmodel_DMDc,Hu,tspan(1:Nt),Hx(:,1));

    % Error
    xpred = zeros(length(Ndelay:Nt+Ndelay-1),2);
    xpred(:,1) = xDMDc(:,3)+repmat(mean(xold(5000:T,1)),[Nt 1]);
    xpred(:,2) = xDMDc(:,4)+repmat(mean(xold(5000:T,2)),[Nt 1]);
    l2err(i) = mean(sqrt(sum((xold(Ndelay:Nt+Ndelay-1,1:2)-xpred).^2,2)));
end

figure,plot(l2err)
[val,idx] = min(l2err)
return
%%
figure,box on,
ccolors = get(gca,'colororder');
plot(tspan,xold(:,1),'-','Color',ccolors(1,:),'LineWidth',1), hold on
plot(tspan,xold(:,2),'-','Color',ccolors(2,:),'LineWidth',1)

plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,3)+repmat(mean(xold(5000:T,1)),[Nt 1]),'.-','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',1)
plot(tspan(Ndelay:Nt+Ndelay-1),xDMDc(:,4)+repmat(mean(xold(5000:T,2)),[Nt 1]),'Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',1)

% plot(tspan(1:Nt),xDMDc(:,1)+repmat(mean(xold(5000:T,1)),[Nt 1]),'.-','Color',ccolors(1,:)-[0 0.2 0.2],'LineWidth',1)
% plot(tspan(1:Nt),xDMDc(:,2)+repmat(mean(xold(5000:T,2)),[Nt 1]),'Color',ccolors(2,:)-[0.1 0.2 0.09],'LineWidth',1)

xlim([0 100])
xlabel('Time')
ylabel('Population size')
legend('Prey (True)','Predator (True)', 'Prey (DMDc)','Predator (DMDc)')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
% print('-depsc2', [figpath,'EX_LOTKA_ModelIdentification_DMDc.eps']);
