

% Parameters True Model
p.SIGMA = 10;             
p.RHO = 28;
p.BETA = 8/3;
Nvar = 3;

% Get training data
ONLY_TRAINING_LENGTH = 1;
InputSignalType = 'sphs';
getTrainingData

% Get validation data
InputSignalType = 'sine2';
getValidationData

% Reference
forcing = @(x,t) (5*sin(30*t)).^3;
[tA,xA] = ode45(@(t,x)LorenzSys(t,x,forcing(x,t),p),tspanV,x(end,1:3),options);   
uA = forcing(0,tspan);

% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
N  = 10;                        % Control / prediction horizon (number of iterations)
Duration = 3;                   % Run for 'Duration' time units
Ton = 0;                        % Time units when control is turned on
getMPCparams
          
x0n=xA(end,:)';                 % Initial condition
% xref1 = [sqrt(p.BETA*(p.RHO-1));sqrt(p.BETA*(p.RHO-1));p.RHO-1];     % critical point
xref1 = [-sqrt(p.BETA*(p.RHO-1));-sqrt(p.BETA*(p.RHO-1));p.RHO-1];
Tcontrol = tA(end);             % Time offset to combine training, prediction and control phase



% Parameters Models
ModelCollection = {'DMDc','DelayDMDc', 'SINDYc', 'NARX'};

clear Results
Results(1:Nmodels) = struct('x',[], 'u', [], 't', [], 'xref', [], 'J', []);

%% Execute MPC for all models
for iM = 1:Nmodels
    disp(['Select model ',num2str(iM) ' of ',num2str(Nmodels)])
    
    this_model = Models(iM);
    
    Ts = this_model.dt; % Sampling time
    Nt = (Duration/Ts)+1;
    
    % Prepare variables
    uopt0    = 0;
    xhat     = x0n;
    uopt     = uopt0.*ones(N,1);
    xHistory = zeros(Nvar,Nt); xHistory(:,1) = xhat;
    uHistory = zeros(1,Nt); uHistory(1)   = uopt(1);
    tHistory = zeros(1,Nt); tHistory(1)   = Tcontrol;
    
    % Parameters Model
    switch select_model
        case 'DelayDMDc'
        case 'DMDc'
        case 'SINDYc'
            this_model.ahat = this_model.Xi(:,1:Nvar);
            
        case 'NARX'
            this_model.Ndelay = 1;
    end
    
    % Start simulation
    fprintf('Simulation started.  It might take a while...\n')
    tic
    try
        for ct = 1:(Duration/Ts)
            % Set references
            xref = xref1;
            
            % NMPC with full-state feedback
            COSTFUN = @(u) ObjectiveFCN_models(u,xhat,N,xref,uopt(1),this_model,diag(Q),R,Ru,select_model);
            uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options);
            
            
            % Integrate system
            xhat = rk4u(@lorenzcontrol_discrete,xhat,uopt(1),Ts/10,10,[],p); %10, 2
            xHistory(:,ct+1) = xhat;
            uHistory(:,ct+1) = uopt(1);
            tHistory(:,ct+1) = ct*Ts+Tcontrol;
             
        end
        fprintf('Simulation finished!\n')
    catch
        fprintf('Simulation ended early!\n')
    end
    %
    tElapsed = toc
    
    Results(iM).eval_time = tElapsed;
    Results(iM).xref = repmat(xref,[1 length(tHistory)]);
    Results(iM).x = xHistory;
    Results(iM).u = uHistory;
    Results(iM).t = tHistory;
    Results(iM).J = evalObjectiveFCN(Results(iM).u,Results(iM).x,Results(iM).xref,diag(Q),R,Ru);
end

%% Unforced
[tU,xU] = ode45(@(t,x)LorenzSys(t,x,0,p),[Tcontrol:Ts:Tcontrol+Duration],xA(end,1:3),options);   % true model