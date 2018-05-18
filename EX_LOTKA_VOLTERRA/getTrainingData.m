%% Execute this file to collect training data for model identification


% Parameters: Model
a = .5;
b = .025;
d = .005;
g = .5;
n = 2;
x0=[60,50];
xref = [g/d;a/b]; % critical point


dt = 0.1;



if ONLY_TRAINING_LENGTH == 1
    % Training
    tspan=[0:dt:200];
    Ntrain = 1001;%ceil(length(tspan)/2);
else
    % Noise
    tspan=[0:dt:400];
    Ntrain = 3001;%ceil(length(tspan)/2);
end

% Time
% t=[dt:dt:100];


options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));

switch InputSignalType
    case 'sine2'
        A = 2;
        forcing = @(x,t) [(A*(sin(1*t)+sin(.1*t))).^2];
        [t,x]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);
        u = forcing(0,tspan);
        
    case 'chirp'
        A = 2;
        forcing = @(x,t) A*chirp(t,[],max(tspan),1.).^2;
        %         forcing = @(x,t) 1*chirp(t,[],100-2*roundn(t,1)+0.01,1.);
        %         forcing = @(x,t) Amp*sawtooth(t,1).*chirp(t,[],max(tspan),1.5); %3.5
        [t,x]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);
        u = forcing(0,tspan);
        
    case 'noise'
        vareps = 0.01; %0.2
        %         forcing = @(x,t) Amp*randn(1,1);
        Diff = @(t,x) [0; vareps];
        SDE = sde(@(t,x) lotkacontrol(t,x,0,a,b,d,g),Diff,'StartState',x0');
        rng(1,'twister')
        %[xdat,t] = simulate(SDE, N, 'DeltaTime', dt);
        [x, t, u] = simByEuler(SDE, length(tspan), 'DeltaTime', dt);
        u = u';
        x = x(1:end-1,:); t = t(1:end-1);
        %         noise = @(t,A)[0;0.0001];
        %         opts = sdeset('RandSeed',19);
        %         [x,u] = sde_euler(lotkacontrol(t,x,0,a,b,d,g),noise,tspan,x0,opts);
        %         t = tspan;
    case 'prbs'
        A = 1; %0.5;
        taulim = [0.1 5];
        states = [-1 1];
        Nswitch = 300;
        forcing = @(x,t) A*prbs(taulim, Nswitch, states, t,0);
        
        [t,x]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);
        
        u = zeros(size(tspan));
        for i = 1:length(tspan)
            u(i) = forcing(0,tspan(i));
        end
        figure,plot(tspan,u)
        
    case 'sphs'
        Pf = 20; %5 % Fundamental period
        K = 8; %16
        A = 5;
        forcing = @(x,t) A*sphs(Pf,K,t).^2;
        [t,x]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);
        u = forcing(0,tspan);
        
    case 'mixed'
        A = [4 5 4 3];
        Pf = 8; % Fundamental period
        K = 16;
        taulim = [0.1 3];
        states = [-1 1];
        Nswitch = 40;
        
        tspan1 = [0:dt:100];
        forcing = @(x,t) SI_Input(t,A,Pf,K,tspan1,taulim, Nswitch, states);
        [t,x] = ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan1,x0,options);
        
        u = zeros(size(tspan1));
        for i = 1:length(tspan1)
            u(i) = forcing(0,tspan1(i));
        end
        figure,plot(tspan1,u)
        
        
        A = 2;
        tspanv = [100:dt:200];
        forcingV = @(x,t) [(A*(sin(1*t)+sin(.1*t))).^2];
        [tv,xv]=ode45(@(t,x) lotkacontrol(t,x,forcingV(x,t),a,b,d,g),tspanv,x(end,:),options);
        uv = forcingV(0,tspanv);
        
        u = [u, uv(2:end)];
        x = [x; xv(2:end,:)];
        t = tspan';
end


%% Split into training and validation data set
%Ntrain = ceil(length(tspan)/2);
xv = x(Ntrain+1:end,:);
x = x(1:Ntrain,:);

uv = u(Ntrain+1:end);
u = u(1:Ntrain);

tv = t(Ntrain+1:end);
t = t(1:Ntrain);

tspanv = tspan(Ntrain+1:end);
tspan = tspan(1:Ntrain);

figure,plot(tspan,u)

%%

figure;
plot(t,x,'LineWidth',1.5)
xlabel('Time')
ylabel('Population size')
legend('Prey','Predator')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
% print('-depsc2', [figpath,'EX_LOTKA_Dynamics.eps']);


T = length(tspan);
N = length(tspan);
% xold = x;
% xmean = mean(xold(5000:T,:));
% xmean = mean(xold);
xmean = xref';

%% Time delay data
if exist('Ndelay','var') == 1
    X   = x - repmat(xmean,[T 1]);
    Hx  = getHankelMatrix_MV(X,max(Ndelay));
    if length(Ndelay) == 1 && Ndelay>1
        Hx = Hx([1:2,2*Ndelay-1:2*Ndelay],:);
    elseif length(Ndelay) == 1 && Ndelay==1
        Hx = Hx(1:2,:);
    end
    Nt  = size(Hx,2);
    
    Hu = getHankelMatrix_MV(u',max(Ndelay));
    if length(Ndelay) == 1 && Ndelay>1
        Hu = Hu([1,Ndelay],:);
    elseif length(Ndelay) == 1 && Ndelay==1
        Hu = Hu;
    end
    % Hu = getHankelMatrix_MV(u',1); Hu = Hu(:,1:Nt);
end