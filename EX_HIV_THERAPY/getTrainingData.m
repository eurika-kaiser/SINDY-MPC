%% Execute this file to collect training data for model identification

%% Parameters: Model
run_HIV_params
rng(3,'twister')
%% Initial Condition
xSTEADY1 = zeros(1,5);
xSTEADY1(2) = ((c2*(lambda1-d*q)-b2*alpha1) - ...
    sqrt((c2*(lambda1-d*q)-b2*alpha1)^2 - 4*alpha1*c2*q*d*b2))/(2*alpha1*c2*q);
xSTEADY1(1) = lambda1/(d+alpha1*xSTEADY1(2));
xSTEADY1(4) = 0;
xSTEADY1(5) = (xSTEADY1(2)*c2*(alpha1*q-a) + b2*alpha1)/(c2*p2*xSTEADY1(2));
xSTEADY1(3) = h*xSTEADY1(5)/(c2*q*xSTEADY1(2));

xSTEADY2 = zeros(1,5);
xSTEADY2(1) = (lambda1*c1)/(d*c1 + b1*alpha1);
xSTEADY2(2) = b1/c1;
xSTEADY2(3) = 0;
xSTEADY2(4) = (alpha1*xSTEADY2(1)-a)/p1;
xSTEADY2(5) = 0;

% KWON
% xSTEADY1 = [lambda1/d, 0, 0, 0, 0];
% xSTEADY2 = zeros(1,5); 
% % A = mu*k1*alpha1/(1-k2*alpha2); % from paper, but wrong
% A = (mu/(1-mu)) / (1-(k2*alpha2)/(k1*alpha1*(1-mu)));
% xSTEADY2(1) = delta*c/(k1*alpha1*(1-mu)); % infected state (steady state)
% xSTEADY2(2) = c*(lambda1-d*xSTEADY2(1))/(xSTEADY2(1)*(k1*alpha1+k2*alpha2*A));
% xSTEADY2(3) = xSTEADY2(2)*A;
% xSTEADY2(4) = alpha1*xSTEADY2(2)/c; 
% xSTEADY2(5) = alpha2*xSTEADY2(3)/c; 
% x0 = [lambda1/d/10, 0, 0, 0.3, 0]; % slightly infected state (close to steady state)

x0 = [10, 0.1, 0.1, 1, 0.1];

n = length(x0);
dt = 1/12;%1/12;%0.01; %0.005;%0.0025;%1/256;%1/12;
xref = xSTEADY1; % healthy steady state


tspan=[0:dt:200]; %200
Ntrain = (length(tspan)-1)/2+1;


options = odeset('RelTol',1e-14,'AbsTol',1e-14*ones(1,n));

switch InputSignalType
    case 'unforced'
        forcing = @(x,t) [1.0];
        [t,x]=ode45(@(t,x) HIVsys_ZURAKOWSKI(t,x,forcing(x,t)),tspan,x0,options);
%         u = forcing(0,tspan)'; 
        u = zeros(length(tspan),1);
        for i = 1:length(tspan)
            u(i,:) = forcing(0,tspan(i));
        end
    case 'sine2'
        A = 1.5;
        Nic = 1; iIC = 1;
%         Nrand = ones(1,5);%
%         Nrand = [1*rand(Nic,9)];
        Nrand = [rand(Nic,1),rand(Nic,1),randn(Nic,1),rand(Nic,1),randn(Nic,1),rand(Nic,1),randn(Nic,1),rand(Nic,1),randn(Nic,1)];
        forcing = @(x,t) [mod((A* (sin(Nrand(iIC,2)*0.7*(t-Nrand(iIC,3))).* ...
            sin(Nrand(iIC,4)*.1*(t-Nrand(iIC,5))).*sin(Nrand(iIC,6)*.2*(t-Nrand(iIC,7))).*sin(Nrand(iIC,8)*.05*(t-Nrand(iIC,9)))) ).^2,1)];
%         forcing = @(x,t) [abs(A* (sin(0.7*t).*sin(.1*t).*sin(.2*t).*sin(.05*t))) ];
%         forcing = @(x,t) [(0.4*(sin(2*pi*0.2*t).*sin(2*pi*0.1*t)))+0.3];
%         forcing = @(x,t) [(0.4*(sin(2*pi*0.2*t).*sin(2*pi*0.05*t)))+0.3];
        [t,x]=ode45(@(t,x) HIVsys_ZURAKOWSKI(t,x,forcing(x,t)),tspan,x0,options);
        u = forcing(0,tspan)';
        
    case 'prbs'
        A = 1; % PAPER
        taulim = [0.2 10];
        states = [0:0.25:1,0,0,0,0];
        Nswitch = 200;
        
%         A = 1;%0.7;%[0.7,0.4]; %0.5;
%         taulim = [0.2 20];%[0.2 10];%[1 5];  %[0.1 0.2]; 
%         states = [0,0,0,0,0:0.25:1,0,0,0,0];%[0,0,0,0.5,0,1,0,1,0,0,0];%[0:0.25:1,0,0,0,0];
        Nswitch = 200;
        forcing = @(x,t) [A(1)*prbs(taulim, Nswitch, states, t,0)];
        
        [t,x]=ode45(@(t,x) HIVsys_ZURAKOWSKI(t,x,forcing(x,t)),tspan,x0,options);
        
        u = zeros(length(tspan),1);
        for i = 1:length(tspan)
            u(i,:) = forcing(0,tspan(i));
        end
        figure,plot(tspan,u)
        size(u)
        
    case 'sphs'
        Pf = 5; %5 % Fundamental period
        K = 3; %16
        A = .1;
        forcing = @(x,t) A*sphs(Pf,K,t).^2;
        [t,x]=ode45(@(t,x) HIVsys_KWON(t,x,forcing(x,t)),tspan,x0,options);
        u = forcing(0,tspan)';
        
end

% x = x(200/dt+1:end,:);
% u = u(200/dt+1:end);
% t = t(200/dt+1:end);
% tspan = t;

%% Split into training and validation data set
%Ntrain = ceil(length(tspan)/2);
% xv = x(Ntrain+1:end,:);
% x = x(1:Ntrain,:);
% 
% uv = u(Ntrain+1:end,:);
% u = u(1:Ntrain,:);
% 
% tv = t(Ntrain+1:end);
% t = t(1:Ntrain);
% 
% tspanv = tspan(Ntrain+1:end);
% tspan = tspan(1:Ntrain);

xv = x;
uv = u;
tv = t;
tspanv = tspan;

figure,plot(tspan,u)

%%
betta_eff = alpha1*(1-eta*u(1));
betta_eff<c1*(c2*b1*(lambda1-q*d)-b2*c1*d)/(b1*(c2*b1*q+b2*c1))
xSTEADY = xSTEADY1;
figure; 
subplot(5,1,1),hold on
% plot(t,xSTEADY(1)*ones(size(t)),'--k')
plot(t,x(:,1),'LineWidth',1.5)
ylabel('x1')
set(gca,'LineWidth',1, 'FontSize',14)
axis tight

subplot(5,1,2),hold on
% plot(t,xSTEADY(2)*ones(size(t)),'--k')
plot(t,x(:,2),'LineWidth',1.5)
ylabel('x2')
set(gca,'LineWidth',1, 'FontSize',14)
axis tight

subplot(5,1,3),hold on
% plot(t,xSTEADY(3)*ones(size(t)),'--k')
plot(t,x(:,3),'LineWidth',1.5)
ylabel('x3')
set(gca,'LineWidth',1, 'FontSize',14)
axis tight

subplot(5,1,4),hold on
% plot(t,xSTEADY(4)*ones(size(t)),'--k')
plot(t,x(:,4),'LineWidth',1.5)
ylabel('x4')
set(gca,'LineWidth',1, 'FontSize',14)
axis tight

subplot(5,1,5),hold on
% plot(t,xSTEADY(5)*ones(size(t)),'--k')
plot(t,x(:,5),'LineWidth',1.5)
ylabel('x5')
set(gca,'LineWidth',1, 'FontSize',14)
axis tight

T = length(tspan);
N = length(tspan);
