%% Execute this file to collect training data for model identification

% Parameters: Model
p.SIGMA = 10;             % True system parameters
p.RHO = 28;
p.BETA = 8/3;
n = 3;
x0=[-8; 8; 27];  % Initial condition
xref1 = [sqrt(p.BETA*(p.RHO-1));sqrt(p.BETA*(p.RHO-1));p.RHO-1];     % critical point
xref2 = [-sqrt(p.BETA*(p.RHO-1));-sqrt(p.BETA*(p.RHO-1));p.RHO-1];

dt = 0.001;  % Time step


if ONLY_TRAINING_LENGTH == 1
    % Without noise-corruption
    tspan=[0:dt:20];
    Ntrain = (length(tspan)-1)/2+1;
else
    % If Noise-corruption
    tspan=[0:dt:40];
    Ntrain = 3*(length(tspan)-1)/4+1;
end


options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));

switch InputSignalType
    case 'sine2'
        A = 2;
        forcing = @(x,t) [(A*(sin(1*t)+sin(.1*t))).^2];
        [t,x]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspan,x0,options);
        u = forcing(0,tspan);
        
    case 'chirp'
        A = 2;
        forcing = @(x,t) A*chirp(t,[],max(tspan),1.).^2;
        [t,x]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspan,x0,options);
        u = forcing(0,tspan);
        
    case 'noise'
        vareps = 0.01; 
        Diff = @(t,x) [0; vareps];
        SDE = sde(@(t,x) LorenzSys(t,x,forcing(x,t),p),Diff,'StartState',x0');
        rng(1,'twister')
        [x, t, u] = simByEuler(SDE, length(tspan), 'DeltaTime', dt);
        u = u';
        x = x(1:end-1,:); t = t(1:end-1);
    case 'prbs'
        A = 1; 
        taulim = [0.1 5];
        states = [-1 1];
        Nswitch = 300;
        forcing = @(x,t) A*prbs(taulim, Nswitch, states, t,0);
        
        [t,x]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspan,x0,options);
        
        u = zeros(size(tspan));
        for i = 1:length(tspan)
            u(i) = forcing(0,tspan(i));
        end
        figure,plot(tspan,u)
        
    case 'sphs'
        Pf = 0.5; % Fundamental period
        K = 8; 
        A = 10;
        forcing = @(x,t) A*sphs(Pf,K,t);
        [t,x]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspan,x0,options);
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
        [t,x] = ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspan1,x0,options);
        
        u = zeros(size(tspan1));
        for i = 1:length(tspan1)
            u(i) = forcing(0,tspan1(i));
        end
        figure,plot(tspan1,u)
        
        
        A = 2;
        tspanv = [100:dt:200];
        forcingV = @(x,t) [(A*(sin(1*t)+sin(.1*t))).^2];
        [tv,xv]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspanv,x(end,:),options);
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

T = length(tspan);
N = length(tspan);

%% Show data
figure,plot(tspan,u)
figure;
plot(t,x,'LineWidth',1.5)
xlabel('Time')
ylabel('Population size')
legend('Prey','Predator')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')


