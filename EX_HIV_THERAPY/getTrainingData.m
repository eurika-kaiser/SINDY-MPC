%% Execute this file to collect training data for model identification

%% Parameters: Model
n = 5;
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
SIM_DURATION = 10;  % Physical time simulation is integrated before integrator is switched, adjust if necessary, depends on length of time series to be computed
dt = 1/12;       % time step of data

run_HIV_params

% Steady-state of model corresponding to progressive infection
xSTEADY1 = zeros(1,5);
xSTEADY1(2) = ((c2*(lambda1-d*q)-b2*alpha1) - ...
    sqrt((c2*(lambda1-d*q)-b2*alpha1)^2 - 4*alpha1*c2*q*d*b2))/(2*alpha1*c2*q);
xSTEADY1(1) = lambda1/(d+alpha1*xSTEADY1(2));
xSTEADY1(4) = 0;
xSTEADY1(5) = (xSTEADY1(2)*c2*(alpha1*q-a) + b2*alpha1)/(c2*p2*xSTEADY1(2));
xSTEADY1(3) = h*xSTEADY1(5)/(c2*q*xSTEADY1(2));

% Steady-state of model corresponding to successful immune response
xSTEADY2 = zeros(1,5);
xSTEADY2(1) = (lambda1*c1)/(d*c1 + b1*alpha1);
xSTEADY2(2) = b1/c1;
xSTEADY2(3) = 0;
xSTEADY2(4) = (alpha1*xSTEADY2(1)-a)/p1;
xSTEADY2(5) = 0;

xref = xSTEADY1; % recovered healthy steady state as reference, might be used in system identification
%% Collect data
if DATA_ENSEMBLE == 1
    
    Nic = 32; %568
    if exist(fullfile(datapath,['DATA_',SystemModel,'_TRAINING-ENSEMBLE_',InputSignalType,'_N',num2str(Nic),'.mat']))==2
        load(fullfile(datapath,['DATA_',SystemModel,'_TRAINING-ENSEMBLE_',InputSignalType,'_N',num2str(Nic),'.mat']))
    else
        tspan=[0:dt:20];
        Ntrain = (length(tspan)-1)/2+1;
        Nt = length(tspan);
        
        if strcmp(ModelName,'SINDYc')
            [x10, x20, x30, x40, x50] = ndgrid([1000],[0,10], [0,10], [0.1,1,10], [0.1,1,10]);
        elseif strcmp(ModelName,'NARX') || strcmp(ModelName,'DMDc')
            [x10, x20, x30, x40, x50] = ndgrid([10,10],[0.1,0.1], [0.1,0.1],[0.1,0.1], [0.1,0.1]);
            
        end
        [N1,N2,N3,N4,N5] = size(x10);
        x0_ensemble = [reshape(x10,[N1*N2*N3*N4*N5,1]), reshape(x20,[N1*N2*N3*N4*N5,1]), reshape(x30,[N1*N2*N3*N4*N5,1]), reshape(x40,[N1*N2*N3*N4*N5,1]), reshape(x50,[N1*N2*N3*N4*N5,1])];
        Nic = size(x0_ensemble,1)
        
        xensemble = zeros(Nt,n,Nic); % Init ensemble data 
        switch InputSignalType
             
            case 'sine2'
                
                rng(1,'twister')
                Nrand = [1*rand(Nic,5)];
                A = 2;
                forcing = @(x,t) [(0.4*(sin(2*pi*0.2*t).*sin(2*pi*0.05*t)))+0.3];
                for iIC = 1:Nic
                    tic
                    try
                        forcing = @(x,t) [(Nrand(iIC,1)*A* (sin(Nrand(iIC,1)*0.7*t).*sin(Nrand(iIC,1)*.1*t).*sin(Nrand(iIC,1)*.2*t).*sin(Nrand(iIC,1)*.05*t)) ).^2]; %[(A*(sin(0.01*t)+sin(.1*t))).^2];
                        [t,x]=ode45(@(t,x) HIVsys_KWON(t,x,forcing(x,t)),tspan,x0_ensemble(iIC,:),options);
                        xensemble(:,:,iIC) = x;
                        for i = 1:length(tspan)
                            u(i,:,iIC) = forcing(0,tspan(i));
                        end
                    catch
                        disp(['ERROR: Simulation failed for i=',num2str(iIC)])
                        xensemble(:,:,iIC) = nan(Nt,n);
                        u(:,:,iIC) = zeros(Nt,length(A));
                    end
                    tend = toc;
                    disp(['Time for ',num2str(iIC), ' of ',num2str(Nic),': ', num2str(tend)])
                end
                figure,plot(t,xensemble(:,:,iIC))
                figure,plot(t,squeeze(xensemble(:,5,:)))
                
            case 'prbs'
                A = 1;
                taulim = [0.2 8];
                states = [0,0:0.25:1,0,0,0]
                Nswitch = 200;
                
                rng(1,'twister');
                seed = randi(Nic,Nic,1); % vary actuation in each trajectory
                forcing = @(x,t,seedval) [A(1)*prbs(taulim, Nswitch, states, t,0,seedval)];
                u = zeros(length(tspan),length(A),Nic);
                
                for iIC = 1:Nic
                    tic
                    try
                        % IF simulation takes longer than SIM_DURATION,
                        % stop simulation and switch solver
                        [t,x,isterminal] = integrateODE(tspan, x0_ensemble(iIC,:), 'ode45', SIM_DURATION,@(x,t) forcing(x,t,seed(iIC)));
                        if isterminal == 1 % try different solver
                            [t,x,isterminal] = integrateODE(tspan, x0_ensemble(iIC,:), 'ode15s', SIM_DURATION,@(x,t) forcing(x,t,seed(iIC)));
                        end
                        xensemble(:,:,iIC) = x;
                        for i = 1:length(tspan)
                            u(i,:,iIC) = forcing(0,tspan(i),seed(iIC));
                        end
                    catch
                        disp(['ERROR: Simulation failed for i=',num2str(iIC)])
                        xensemble(:,:,iIC) = nan(Nt,n);
                        u(:,:,iIC) = zeros(Nt,length(A));
                    end
                    tend = toc;
                    disp(['Time for ',num2str(iIC), ' of ',num2str(Nic),': ', num2str(tend)])
                end
                
        end
        
        
        %% Clean up data
        % from unstable simulations or which failed
        IXnaninf = zeros(Nic,1);
        for i = 1:Nic
            IXnan = isnan(squeeze((xensemble(:,:,i))));
            IXinf = isinf(squeeze(abs(xensemble(:,:,i))));
            if any(IXnan(:)==1) || any(IXinf(:)==1)
                IXnaninf(i) = 1;
            end
        end
        xensemble(:,:,logical(IXnaninf)) = [];
        u(:,:,logical(IXnaninf)) = [];
        Nic = size(xensemble,3);
        
        %% Split into training and validation data set
        Ntrain = ceil(length(tspan)/2);
        xv = xensemble(Ntrain+1:end,:,:);
        x = xensemble(1:Ntrain,:,:);
        
        uv = u(Ntrain+1:end,:,:);
        u = u(1:Ntrain,:,:);
        
        tv = t(Ntrain+1:end);
        t = t(1:Ntrain);
        
        tspanv = tspan(Ntrain+1:end);
        tspan = tspan(1:Ntrain);
       
        
        %% Show data
        
        figure;
        subplot(6,1,1)
        plot(tspan,squeeze(x(:,1,:)),'LineWidth',1.5)
        ylabel('x1')
        set(gca,'LineWidth',1, 'FontSize',14)
        
        subplot(6,1,2)
        plot(tspan,squeeze(x(:,2,:)),'LineWidth',1.5)
        ylabel('x2')
        set(gca,'LineWidth',1, 'FontSize',14)
        
        subplot(6,1,3)
        plot(tspan,squeeze(x(:,3,:)),'LineWidth',1.5)
        ylabel('x3')
        set(gca,'LineWidth',1, 'FontSize',14)
        
        subplot(6,1,4)
        plot(tspan,squeeze(x(:,4,:)),'LineWidth',1.5)
        ylabel('x4')
        set(gca,'LineWidth',1, 'FontSize',14)
        
        subplot(6,1,5)
        plot(tspan,squeeze(x(:,5,:)),'LineWidth',1.5)
        ylabel('x5')
        set(gca,'LineWidth',1, 'FontSize',14)
        
        subplot(6,1,6)
        plot(tspan,squeeze(u(:,1,:)),'LineWidth',1.5)
        ylabel('u')
        set(gca,'LineWidth',1, 'FontSize',14)
        %%
        T = length(tspan);
        N = length(tspan);
        save(fullfile(datapath,['DATA_',SystemModel,'_TRAINING-ENSEMBLE_',InputSignalType,'_N',num2str(Nic),'.mat']))
    end
elseif DATA_ENSEMBLE == 0
    rng(3,'twister')
    
    % Initial Condition
    x0 = [10, 0.1, 0.1, 1, 0.1];
    tspan=[0:dt:400];
    Ntrain = (length(tspan)-1)/2+1;
    
   
    
    switch InputSignalType
        case 'unforced'
            forcing = @(x,t) [0];
            [t,x]=ode45(@(t,x) HIVsys_ZURAKOWSKI(t,x,forcing(x,t)),tspan,x0,options);
            u = zeros(length(tspan),1);
            for i = 1:length(tspan)
                u(i,:) = forcing(0,tspan(i));
            end
        case 'sine2'
            A = 1.5;
            Nic = 1; iIC = 1;
            Nrand = [rand(Nic,1),rand(Nic,1),randn(Nic,1),rand(Nic,1),randn(Nic,1),rand(Nic,1),randn(Nic,1),rand(Nic,1),randn(Nic,1)];
            forcing = @(x,t) [mod((A* (sin(Nrand(iIC,2)*0.7*(t-Nrand(iIC,3))).* ...
                sin(Nrand(iIC,4)*.1*(t-Nrand(iIC,5))).*sin(Nrand(iIC,6)*.2*(t-Nrand(iIC,7))).*sin(Nrand(iIC,8)*.05*(t-Nrand(iIC,9)))) ).^2,1)];
            [t,x]=ode45(@(t,x) HIVsys_ZURAKOWSKI(t,x,forcing(x,t)),tspan,x0,options);
            u = forcing(0,tspan)';
            
        case 'prbs'
            A = 1; % PAPER
            taulim = [0.2 10];
            states = [0:0.25:1,0,0,0,0];
            Nswitch = 200;
            forcing = @(x,t) [A(1)*prbs(taulim, Nswitch, states, t,0)];
            
            [t,x]=ode45(@(t,x) HIVsys_ZURAKOWSKI(t,x,forcing(x,t)),tspan,x0,options);
            
            u = zeros(length(tspan),1);
            for i = 1:length(tspan)
                u(i,:) = forcing(0,tspan(i));
            end
            figure,plot(tspan,u)
            size(u)
            
    end
    
    
    %% Split into training and validation data set
    Ntrain = ceil(length(tspan)/2);
    xv = x(Ntrain+1:end,:);
    x = x(1:Ntrain,:);
    
    uv = u(Ntrain+1:end,:);
    u = u(1:Ntrain,:);
    
    tv = t(Ntrain+1:end);
    t = t(1:Ntrain);
    
    tspanv = tspan(Ntrain+1:end);
    tspan = tspan(1:Ntrain);
    
    %% Show results
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
end