%% Execute this file to collect validation data for model identification

x0 = x(end,:);
tspanV =[tspan(end):dt:20];

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));

switch InputSignalType
    case 'sine2'
        A = 2;
        forcing = @(x,t) [(A*(sin(1*t)+sin(.1*t))).^2]; 
        [t_valid,x_valid]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspanV,x0,options);
        u_valid = forcing(0,tspanV);
        
    case 'chirp'
        A = 10;
        forcing = @(x,t) A*chirp(t,[],max(tspanV),1.).^2;
        [t_valid,x_valid]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspanV,x0,options);
        u_valid = forcing(0,tspanV);
        
    case 'noise'
        vareps = 0.01; 
        Diff = @(t,x) [0; vareps];
        SDE = sde(@(t,x) LorenzSys(t,x,forcing(x,t),p),Diff,'StartState',x0');
        rng(1,'twister')
        [x_valid, t_valid, u_valid] = simByEuler(SDE, length(tspanV), 'DeltaTime', dt);
        u_valid = u_valid';
        x_valid = x_valid(1:end-1,:); t_valid = t_valid(1:end-1);
    case 'prbs'
        A = 1; 
        taulim = [0.1 5];
        states = [-1 1];
        Nswitch = 300;
        forcing = @(x,t) A*prbs(taulim, Nswitch, states, t,0);
        
        [t_valid,x_valid]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspanV,x0,options);
        
        u_valid = zeros(size(tspanV));
        for i = 1:length(tspanV)
            u_valid(i) = forcing(0,tspanV(i));
        end
        figure,plot(tspanV,u_valid)
        
    case 'sphs'
        Pf = 1; % Fundamental period
        K = 8; 
        A = 10;
        forcing = @(x,t) A*sphs(Pf,K,t);
        [t_valid,x_valid]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspanV,x0,options);
        u_valid = forcing(0,tspanV);
        
    case 'mixed'
        A = [4 5 4 3];
        Pf = 8; % Fundamental period
        K = 16;
        taulim = [0.1 3];
        states = [-1 1];
        Nswitch = 40;
        
        tspan1 = [0:dt:100];
        forcing = @(x,t) SI_Input(t,A,Pf,K,tspan1,taulim, Nswitch, states);
        [t_valid,x_valid] = ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspan1,x0,options);
        
        u_valid = zeros(size(tspan1));
        for i = 1:length(tspan1)
            u_valid(i) = forcing(0,tspan1(i));
        end
        figure,plot(tspan1,u_valid)
        
        
        A = 2;
        tspanv = [100:dt:200];
        forcingV = @(x,t) [(A*(sin(1*t)+sin(.1*t))).^2]; 
        [tv,xv]=ode45(@(t,x) LorenzSys(t,x,forcing(x,t),p),tspanv,x_valid(end,:),options);
        uv = forcingV(0,tspanv);
        
        u_valid = [u_valid, uv(2:end)];
        x_valid = [x_valid; xv(2:end,:)];
        t_valid = tspanV';
end

T = length(tspanV);
N = length(tspanV);

%% Show data
figure;
plot(t_valid,x_valid,'LineWidth',1.5)
xlabel('Time')
ylabel('xi')
legend('x','y','z')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')



