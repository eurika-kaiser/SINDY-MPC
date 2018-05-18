%% Execute this file to collect training data for model identification

% Parameters: Model
x0 = [0.3 0.1 -0.5]; % Initial condition
n = 3;
dt = 0.01;

if ONLY_TRAINING_LENGTH == 1
    % Training
    tspan=[0:dt:40]; %60
    Ntrain = (length(tspan)-1)/2+1;
else
    % Noise
    tspan=[0:dt:40];
    Ntrain = 3*(length(tspan)-1)/4+1;
end
Nt = length(tspan);

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));

% Collect data
% [Xgrid,Ygrid,Zgrid] = meshgrid([-0.3:0.1:0.3],[-0.1:0.01:0.1],[-0.5,0.05:0.5]);
% [Xgrid,Ygrid,Zgrid] = meshgrid([-0.5:0.05:0.5],[-0.3:0.05:0.3],[-0.7,0.05:0.7]);
[Xgrid,Ygrid,Zgrid] = meshgrid([-0.5:0.1:0.5],[-0.3:0.1:0.3],[-0.8:0.1:0.8]);
Nx = size(Xgrid,1); Ny = size(Ygrid,2); Nz = size(Zgrid,3);
Nxyz = Nx*Ny*Nz;
y = zeros(Nx*Ny*Nz,Nt,3);
y(:,1,1) = reshape(Xgrid,[Nxyz 1]);
y(:,1,2) = reshape(Ygrid,[Nxyz 1]);
y(:,1,3) = reshape(Zgrid,[Nxyz 1]);

switch InputSignalType
    case 'sine2'
        A = .05;
        forcing = @(x,t) [(A* (sin(0.7*t).*sin(.1*t).*sin(.2*t).*sin(.05*t))) ]; %[(A*(sin(0.01*t)+sin(.1*t))).^2];
        u = forcing(0,tspan);
        
        for k=1:size(y,1)
            [t,tmp] = ode45(@(t,x) F8Sys(t,x,forcing(x,t)),tspan,y(k,1,:),options);
            y(k,:,:) = tmp;
        end
end

%%

figure; box on, hold on
hold on
for k=1:size(y,1)
    plot3(squeeze(y(k,:,1)),squeeze(y(k,:,2)),squeeze(y(k,:,3)),'-k')
end
xlabel('x1:angle of attack'), ylabel('x2:pitch angle'), ylabel('x3:pitch rate')
drawnow
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 220 200])
set(gcf,'PaperPositionMode','auto')



%% Validation data

if ONLY_TRAINING_LENGTH == 1
    % Training
    tspanv=[0:dt:200];
    Ntrain = (length(tspanv)-1)/2+1;
else
    % Noise
    tspanv=[0:dt:400];
    Ntrain = 3*(length(tspanv)-1)/4+1;
end


options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));



switch InputSignalType
    case 'sine2'
        A = .05;
        forcing = @(x,t) [(A* (sin(0.7*t).*sin(.1*t).*sin(.2*t).*sin(.05*t)) )]; %[(A*(sin(0.01*t)+sin(.1*t))).^2];
        [tv,xv]=ode45(@(t,x) F8Sys(t,x,forcing(x,t)),tspanv,x0,options);
        uv = forcing(0,tspanv);    
end

%% Split into training and validation data set
%Ntrain = ceil(length(tspan)/2);
x = y;
u = u;
t = t;

tspanv = tv;
tspan = tspan;

figure,plot(tspan,u)

T = length(tspan);
N = length(tspan);
