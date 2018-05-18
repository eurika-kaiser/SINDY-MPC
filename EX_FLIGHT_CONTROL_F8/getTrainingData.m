%% Execute this file to collect training data for model identification

x0  = [0.3 0.1 -0.5];   % Initial condition
n   = length(x0);       % Number of variables
dt  = 0.01;             % Time step
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));

if ENSEMBLE_DATA == 0
    tspan =[0:dt:200];
    Ntrain = (length(tspan)-1)/2+1;
    switch InputSignalType
        case 'sine2'
            A = .05;
            forcing = @(x,t) [(A* (sin(0.7*t).*sin(.1*t).*sin(.2*t).*sin(.05*t)) )]; %[(A*(sin(0.01*t)+sin(.1*t))).^2];
            [t,x]=ode45(@(t,x) F8Sys(t,x,forcing(x,t)),tspan,x0,options);
            u = forcing(0,tspan);
        case 'sine3'
            forcing = @(x,t) (.5*sin(5*t).*sin(.5*t)+0.1).^3;
            [t,x]=ode45(@(t,x) F8Sys(t,x,forcing(x,t)),tspan,x0,options);
            u = forcing(0,tspan);
        case 'chirp'
            A = .5;
            forcing = @(x,t) A*chirp(t,[],max(tspan),0.1).^2;
            [t,x]=ode45(@(t,x) F8Sys(t,x,forcing(x,t)),tspan,x0,options);
            u = forcing(0,tspan);
        case 'noise'
            vareps = 0.01;
            Diff = @(t,x) [vareps; 0; 0];
            SDE = sde(@(t,x) F8Sys(t,x,0),Diff,'StartState',x0');
            rng(1,'twister')
            [x, t, u] = simByEuler(SDE, length(tspan), 'DeltaTime', dt);
            u = u';
            x = x(1:end-1,:); t = t(1:end-1);
        case 'prbs'
            A = 0.05236; 
            taulim = [1 8];
            states = [-0.5:0.25:0.5];
            Nswitch = 4000;
            forcing = @(x,t) A*prbs(taulim, Nswitch, states, t,0);
            
            [t,x]=ode45(@(t,x) F8Sys(t,x,forcing(x,t)),tspan,x0,options);
            
            u = zeros(size(tspan));
            for i = 1:length(tspan)
                u(i) = forcing(0,tspan(i));
            end
            figure,plot(tspan,u)
            
        case 'sphs'
            Pf = 10;% Fundamental period
            K = 2; 
            A = 0.1;
            forcing = @(x,t) A*sphs(Pf,K,t);
            [t,x]=ode45(@(t,x) F8Sys(t,x,forcing(x,t)),tspan,x0,options);
            u = forcing(0,tspan);
    end
    
    %% Split into training and validation data set
    figure,plot(tspan,u)
    
    figure;
    plot(t,x,'LineWidth',1.5)
    xlabel('Time')
    ylabel('xi')
    legend('angle of attack','pitch angle', 'pitch rate')
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
    
    xv = x(Ntrain+1:end,:);
    x = x(1:Ntrain,:);
    
    uv = u(Ntrain+1:end);
    u = u(1:Ntrain);
    
    tv = t(Ntrain+1:end);
    t = t(1:Ntrain);
    
    tspanv = tspan(Ntrain+1:end);
    tspan = tspan(1:Ntrain);
    
elseif ENSEMBLE_DATA == 1
    tspan =[0:dt:20];
    Ntrain = (length(tspan)-1)/2+1;
    [x10, x20, x30] = ndgrid([-0.4:0.1:0.4],[-4:1:4], [-4:0.5:4]);
    [N1,N2,N3,N4,N5] = size(x10);
    x0_ensemble = [reshape(x10,[N1*N2*N3*N4*N5,1]), reshape(x20,[N1*N2*N3*N4*N5,1]), reshape(x30,[N1*N2*N3*N4*N5,1])];
    Nic = size(x0_ensemble,1);
    Nt = length(tspan);
    xensemble = zeros(Nt,n,Nic);
    u = zeros(Nt,1,Nic);
    switch InputSignalType
        case 'sine2'
            A = .05;
            rng(1,'twister')
            Nrand = [2*rand(Nic,5)];
            
             
        case 'sine3'
            A = .05;
            rng(1,'twister')
            Nrand = [2*rand(Nic,3)];
            
        case 'sphs'
            rng(1,'twister')
            Pf = 1;  % Fundamental period
            K = 5; 
            A = .1;
            Nrand = [rand(Nic,1),2*rand(Nic,2)];
            
        case 'prbs'
            A = 0.05;
            taulim = [1 8];
            states = [-0.3:0.15:0.5,0,0,0,0];
            Nswitch = 4000;
            forcing = @(x,t) A*prbs(taulim, Nswitch, states, t,0);

    end
    
    for iIC = 1:Nic
        tic
        try
            switch InputSignalType
                case 'sphs'
                    forcing = @(x,t) Nrand(iIC,1)*A*sphs(Nrand(iIC,2)*Pf,Nrand(iIC,3)*K,t);
                case 'sine2'
                    forcing = @(x,t) [(Nrand(iIC,5)*A* (sin(Nrand(iIC,1)*0.7*t).*sin(Nrand(iIC,2)*.1*t).*sin(Nrand(iIC,3)*.2*t).*sin(Nrand(iIC,4)*.05*t)) )]; %[(A*(sin(0.01*t)+sin(.1*t))).^2];
                case 'sine3'
                    forcing = @(x,t) (Nrand(iIC,1)*.5*sin(Nrand(iIC,2)*5*t).*sin(Nrand(iIC,3)*.5*t)+0.1).^3;  
            end
            [t,x]=ode45(@(t,x) F8Sys(t,x,forcing(x,t)),tspan,x0_ensemble(iIC,:),options);
            xensemble(:,:,iIC) = x;
            for i = 1:length(tspan)
                u(i,:,iIC) = forcing(0,tspan(i));
            end
        catch
            disp(['ERROR: Simulation failed for i=',num2str(iIC)])
            xensemble(:,:,iIC) = nan(Nt,n);
            u(:,:,iIC) = zeros(Nt,length(A));
        end
        if all(isnan(xensemble(:,1,iIC))==0)
            disp(['DONE: Simulation done for i=',num2str(iIC)])
        end
        tend = toc;
        disp(['Time for ',num2str(iIC), ' of ',num2str(Nic),': ', num2str(tend)])
    end
    
    
    
    %% Split into training and validation data set
    ICinvalid = isnan(squeeze(xensemble(1,1,:)));
    xensemble(:,:,ICinvalid==1) = [];
    u(:,:,ICinvalid==1) = []; 
    
    xv = xensemble(Ntrain+1:end,:,:);
    x = xensemble(1:Ntrain,:,:);
    
    uv = u(Ntrain+1:end,:,:);
    u = u(1:Ntrain,:,:);
    
    tv = tspan(Ntrain+1:end);
    t = tspan(1:Ntrain);
    
    tspanv = tspan(Ntrain+1:end);
    tspan = tspan(1:Ntrain);
    
    Nic = size(xensemble,3);
    figure; hold on,box on
    for i = 1:Nic
        plot(tspan,x(:,1,i),'LineWidth',1.5)
    end
    xlabel('Time')
    ylabel('xi')
    legend('angle of attack','pitch angle', 'pitch rate')
    set(gca,'LineWidth',1, 'FontSize',14)
    set(gcf,'Position',[100 100 300 200])
    set(gcf,'PaperPositionMode','auto')
end

