%% Execute this file to collect training data for model identification

x0 = x(end,:);
tspanV =[tspan(end):dt:200];

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));

switch InputSignalType
    case 'sine2'
        A = 2;
        forcing = @(x,t) [(A*(sin(1*t)+sin(.1*t))).^2]; 
        [t_valid,x_valid]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspanV,x0,options);
        u_valid = forcing(0,tspanV);
        
    case 'chirp'
        A = 2;
        forcing = @(x,t) A*chirp(t,[],max(tspanV),1.).^2;
        [t_valid,x_valid]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspanV,x0,options);
        u_valid = forcing(0,tspanV);
        
    case 'noise'
        vareps = 0.01; 
        Diff = @(t,x) [0; vareps];
        SDE = sde(@(t,x) lotkacontrol(t,x,0,a,b,d,g),Diff,'StartState',x0');
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
        
        [t_valid,x_valid]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspanV,x0,options);
        
        u_valid = zeros(size(tspanV));
        for i = 1:length(tspanV)
            u_valid(i) = forcing(0,tspanV(i));
        end
        figure,plot(tspanV,u_valid)
        
    case 'sphs'
        Pf = 20;  % Fundamental period
        K = 8; 
        A = 5;
        forcing = @(x,t) A*sphs(Pf,K,t).^2;
        [t_valid,x_valid]=ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspanV,x0,options);
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
        [t_valid,x_valid] = ode45(@(t,x) lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan1,x0,options);
        
        u_valid = zeros(size(tspan1));
        for i = 1:length(tspan1)
            u_valid(i) = forcing(0,tspan1(i));
        end
        figure,plot(tspan1,u_valid)
        
        
        A = 2;
        tspanv = [100:dt:200];
        forcingV = @(x,t) [(A*(sin(1*t)+sin(.1*t))).^2]; 
        [tv,xv]=ode45(@(t,x) lotkacontrol(t,x,forcingV(x,t),a,b,d,g),tspanv,x_valid(end,:),options);
        uv = forcingV(0,tspanv);
        
        u_valid = [u_valid, uv(2:end)];
        x_valid = [x_valid; xv(2:end,:)];
        t_valid = tspanV';
end


%%

figure;
plot(t_valid,x_valid,'LineWidth',1.5)
xlabel('Time')
ylabel('Population size')
legend('Prey','Predator')
set(gca,'LineWidth',1, 'FontSize',14)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')

T = length(tspanV);
N = length(tspanV);

%% Time delay data
if exist('Ndelay','var') == 1
    X   = x_valid - repmat(xmean,[T 1]);
    Hx_valid  = getHankelMatrix_MV(X,max(Ndelay));
    if length(Ndelay) == 1 && Ndelay>1
        Hx_valid = Hx_valid([1:2,2*Ndelay-1:2*Ndelay],:);
    elseif length(Ndelay) == 1 && Ndelay==1
        Hx_valid = Hx_valid(1:2,:);
    end
    Nt  = size(Hx_valid,2);
    
    Hu_valid = getHankelMatrix_MV(u_valid',max(Ndelay));
    if length(Ndelay) == 1 && Ndelay>1
        Hu_valid = Hu_valid([1,Ndelay],:);
    elseif length(Ndelay) == 1 && Ndelay==1
        Hu_valid = Hu_valid;
    end
end