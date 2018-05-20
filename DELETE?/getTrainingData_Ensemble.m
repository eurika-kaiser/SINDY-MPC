%% Execute this file to collect training data for model identification
% 
% if exist(fullfile(datapath,['DATA_',SystemModel,'_TRAINING-ENSEMBLE_',InputSignalType,'.mat']))==2
%     load(fullfile(datapath,['DATA_',SystemModel,'_TRAINING-ENSEMBLE_',InputSignalType,'.mat']))
SIM_DURATION = 10;
Nic = 32;%568
if exist(fullfile(datapath,['DATA_',SystemModel,'_TRAINING-ENSEMBLE_',InputSignalType,'_N',num2str(Nic),'.mat']))==2
    load(fullfile(datapath,['DATA_',SystemModel,'_TRAINING-ENSEMBLE_',InputSignalType,'_N',num2str(Nic),'.mat']))
else
    n = 5;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
    %% Parameters: Model
    run_HIV_params
    
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
    % %%%% Initial Condition
    % xSTEADY1 = [lambda1/d, 0, 0, 0, 0];
    %
    % xSTEADY2 = zeros(1,5);
    % % A = mu*k1*alpha1/(1-k2*alpha2); % from paper, but wrong
    % A = (mu/(1-mu)) / (1-(k2*alpha2)/(k1*alpha1*(1-mu)));
    % xSTEADY2(1) = delta*c/(k1*alpha1*(1-mu)); % infected state (steady state)
    % xSTEADY2(2) = c*(lambda1-d*xSTEADY2(1))/(xSTEADY2(1)*(k1*alpha1+k2*alpha2*A));
    % xSTEADY2(3) = xSTEADY2(2)*A;
    % xSTEADY2(4) = alpha1*xSTEADY2(2)/c;
    % xSTEADY2(5) = alpha2*xSTEADY2(3)/c;
    
    % x0 = [lambda1/d, 0, 0, 0.3, 0]; % slightly infected state (close to steady state)
    
    
    % [x10, x20, x30, x40, x50] = ndgrid([1000],[0,1,10], [0,10], [0.1,1,10], [0.1,1,10]); % 972
    if strcmp(ModelName,'SINDYc')
        [x10, x20, x30, x40, x50] = ndgrid([1000],[0,10], [0,10], [0.1,1,10], [0.1,1,10]);
    elseif strcmp(ModelName,'NARX') || strcmp(ModelName,'DMDc')
        %     [x10, x20, x30, x40, x50] = ndgrid([1000],[10], [10], [1,10], [0.1,1,10]);
        %     [x10, x20, x30, x40, x50] = ndgrid([1000],[0,10], [0,10], [0.1,1,10], [0.1,1,10]);
        %     [x10, x20, x30, x40, x50] = ndgrid([1,10,100,1000],[0,1,10], [0,1,10], [0.1,1,10], [0.1,1,10]);
        %     [x10, x20, x30, x40, x50] = ndgrid([100,500,1000],[0,0.1,1,10], [0,0.1,1,10], [0.1,1,10], [0.1,1,10]);
        %     [x10, x20, x30, x40, x50] = ndgrid([10,50,100],[0,.01,10],[0,.01,10], [0,.01,10], [0,.01,10]); %prbs
        %       [x10, x20, x30, x40, x50] = ndgrid([1,10],[0.1,1,10], [0.1,10], [0.1,1,10], [0.1,1,10]);
%             [x10, x20, x30, x40, x50] = ndgrid([10,10],[0.1,0.1], [0.1,0.1],[0.1,0.1], [0.1,0.1]);  % in Paper
            [x10, x20, x30, x40, x50] = ndgrid([5,10,15],[0.05,0.1,0.2], [0.05,0.1,0.2],[0.05,0.1,0.2], [0.05,0.1,0.2]);  % in Paper
%         [x10, x20, x30, x40, x50] = ndgrid([5,10,800,1000],[0.01,0.1], [0.1,1,10,100,1000],[0.1], [0.1,1,8,10]);
%         [x10, x20, x30, x40, x50] =
%         ndgrid([5,10,50,100,200,400,600,800],[0.01,0.1],
%         [0.1,1,10,100,200,400,600,800,1000],[0.1], [0.1,1,8,10]); %568
    end
    [N1,N2,N3,N4,N5] = size(x10);
    x0_ensemble = [reshape(x10,[N1*N2*N3*N4*N5,1]), reshape(x20,[N1*N2*N3*N4*N5,1]), reshape(x30,[N1*N2*N3*N4*N5,1]), reshape(x40,[N1*N2*N3*N4*N5,1]), reshape(x50,[N1*N2*N3*N4*N5,1])];
    Nic = size(x0_ensemble,1),
    % Nic = 1
    % Nic = 30
    
    
    dt = 1/12;%0.01;
    xref = xSTEADY1; % recovered healthy steady state
    
    if ONLY_TRAINING_LENGTH == 1
        % Training
        tspan=[0:dt:20];%200
        Ntrain = (length(tspan)-1)/2+1;
    else
        % Noise
        tspan=[0:dt:20];
        Ntrain = 3*(length(tspan)-1)/4+1;
    end
    
    Nt = length(tspan);
    xensemble = zeros(Nt,n,Nic);
    
    
    switch InputSignalType
        case 'unforced'
            A = 0;
            forcing = @(x,t) [A];
            u = zeros(length(tspan),length(A),Nic);
            for iIC = 1:Nic
                tic
                try
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
            
            
            
            
        case 'sine2'
            A = 2;
            rng(1,'twister')
            Nrand = [1*rand(Nic,5)];
            
            
            %         forcing = @(x,t) [(0.4*(sin(2*pi*0.2*t).*sin(2*pi*0.05*t)))+0.3];
            u = zeros(length(tspan),length(A),Nic);
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
            A = 1;%0.7;%[0.7,0.4]; %0.5;
            taulim = [0.2 8];%[0.2 10];%[1 5];  %[0.1 0.2];
            states = [0,0:0.25:1,0,0,0];%[0,0,0,0.5,0,1,0,1,0,0,0];%[0:0.25:1,0,0,0,0];
            Nswitch = 200;
            
            %         forcing = @(x,t) [A(1)*prbs(taulim, Nswitch, states, t,0)];
            
            rng(1,'twister');
            seed = randi(Nic,Nic,1);
            forcing = @(x,t,seedval) [A(1)*prbs(taulim, Nswitch, states, t,0,seedval)];
            u = zeros(length(tspan),length(A),Nic);
            
            for iIC = 1:Nic
                tic
                try
                    %                 [t,x]=ode45(@(t,x) HIVsys_KWON(t,x,forcing(x,t)),tspan,x0_ensemble(iIC,:),options);
                    %                 [t,x]=ode45(@(t,x) HIVsys_ZURAKOWSKI(t,x,forcing(x,t,seed(iIC))),tspan,x0_ensemble(iIC,:),options);
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
            
            
            %         if exist(['DataEnsemble_N',num2str(Nic),'.mat'])==2
            % %             load('DataEnsemble.mat')
            %             load(['DataEnsemble_N',num2str(Nic),'.mat'])
            %         else
            % %             A = [0.8,0.4]; %0.5;
            %             A = [0.8];
            % %             taulim = [2 10];  %[0.1 0.2];
            %             taulim = [2 10];
            %             states = [0:0.25:1];
            %             Nswitch = 50;
            %
            %             rng(0,'twister');
            %             seed = randi(Nic,Nic,1);
            %             forcing = @(x,t,seedval) [A(1)*prbs(taulim, Nswitch, states, t,0,seedval),A(2)*prbs(taulim, Nswitch, states, t,0,seedval)]; %% ADD SEED!@@@!!!!
            %
            %             u = zeros(length(tspan),length(A),Nic);
            %             for iIC = 1:Nic %1140:Nic
            %                 tic
            %                 try
            %                     [t,x]=ode45(@(t,x) HIVsys(t,x,forcing(x,t,seed(iIC))),tspan,x0_ensemble(iIC,:),options);
            %                     xensemble(:,:,iIC) = x;
            %                     for i = 1:length(tspan)
            %                         u(i,:,iIC) = forcing(0,tspan(i),seed(iIC));
            %                     end
            %                 catch
            %                    disp(['ERROR: Simulation failed for i=',num2str(iIC)])
            %                    xensemble(:,:,iIC) = nan(Nt,n);
            %                    u(:,:,iIC) = zeros(Nt,2);
            %                 end
            %                 tend = toc;
            %                 disp(['Time for ',num2str(iIC), ' of ',num2str(Nic),': ', num2str(tend)])
            %             end
            %             figure,plot(t,xensemble(:,:,iIC))
            %             figure,plot(t,squeeze(xensemble(:,6,:)))
            %
            % %             u = zeros(length(tspan),2);
            % %             for i = 1:length(tspan)
            % %                 u(i,:) = forcing(0,tspan(i));
            % %             end
            %
            %             save(['DataEnsemble_N',num2str(Nic),'.mat'])
            %
            %         end
        case 'sphs'
            Pf = 10; %5 % Fundamental period
            K = 3; %16
            A = .2;
            forcing = @(x,t) A*sphs(Pf,K,t).^2;
            u = zeros(length(tspan),length(A),Nic);
            for iIC = 1:Nic
                tic
                try
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
    end
    
    % x = x(200/dt+1:end,:);
    % u = u(200/dt+1:end);
    % t = t(200/dt+1:end);
    % tspan = t;
    
    %%
    % Nic = iIC-1;
    % xensemble = xensemble(:,:,1:Nic);
    % u = u(:,:,1:Nic);
    
    %% Clean uo data
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
    %Ntrain = ceil(length(tspan)/2);
    % xv = xensemble(Ntrain+1:end,:,:);
    % x = xensemble(1:Ntrain,:,:);
    %
    % uv = u(Ntrain+1:end,:,:);
    % u = u(1:Ntrain,:,:);
    %
    % tv = t(Ntrain+1:end);
    % t = t(1:Ntrain);
    %
    % tspanv = tspan(Ntrain+1:end);
    % tspan = tspan(1:Ntrain);
    
    x = xensemble;
    xv = xensemble;
    
    uv = u;
    
    tv = tspan;
    
    tspanv = tspan;
    
    % figure,plot(tspan,u)
    
    %%
    
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


