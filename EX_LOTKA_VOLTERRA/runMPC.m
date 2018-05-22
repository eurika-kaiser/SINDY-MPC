% Parameters: Model
a = .5;
b = .025;
d = .005;
g = .5;
n = 2;
x0=[60,50];
xref = [g/d;a/b]; % critical point
dt = Models(1).dt;

Pf = 20;  % Fundamental period
K = 8; 
A = 5;
forcing = @(x,t) A*sphs(Pf,K,t).^2;





options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,2));
tspan = [0:dt:200];
[tA,xA] = ode45(@(t,x)lotkacontrol(t,x,forcing(x,t),a,b,d,g),tspan,x0,options);
uA = forcing(0,tspan);

% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
N  = 10;                        % Control / prediction horizon (number of iterations)
Duration = 20;                  % Run for 'Duration' time units
Ton = 0;                        % Time units when control is turned on
Nvar = 2;
Q = [1 1];                      % State weights
R = 0.5;                        % du weights
Ru = 0.5;                       % u weights
B = [0; 1];
C = eye(Nvar);
D = 0;
LB = -20*ones(N,1);             % Lower bound of control input
UB = 20*ones(N,1);              % Upper bound of control input      
x0n=xA(end,:)';                 % Initial condition
xref1 = [g/d;a/b];              % critical point
Tcontrol = tA(end);             % Time offset to combine training, prediction and control phase

% Parameters Models
ModelCollection = {'DMDc','DelayDMDc', 'SINDYc', 'NARX'};
%Nmodels = length(ModelCollection);

% Parameters True Model
p.a = a;
p.b = b;
p.d = d;
p.g = g;

clear Results
Results(1:Nmodels) = struct('x',[], 'u', [], 't', [], 'xref', [], 'J', []);

%%
for iM = 1:Nmodels
    disp(['Select model ',num2str(iM) ' of ',num2str(Nmodels)])
    
    this_model = Models(iM);
    
    Ts = this_model.dt; % Sampling time
    Nt = (Duration/Ts)+1;
    
    % Prepare variables
    uopt0    = 0;
    xhat     = x0n;
    uopt     = uopt0.*ones(N,1);
    xHistory = zeros(2,Nt); xHistory(:,1) = xhat;
    uHistory = zeros(1,Nt); uHistory(1)   = uopt(1);
    tHistory = zeros(1,Nt); tHistory(1)   = Tcontrol;
    
    % Parameters Model
    switch select_model
        case 'DelayDMDc'
            pest.dt = Models.DelayDMDc.dt;
            pest.sys = Models.DelayDMDc.sys;
            pest.xmean = xmean';
            pest.udelay = zeros(1,N);
            pest.xdelay = [x(end-Ndelay+1,1:2)]';
            pest.Nxlim = size(x,1);
            Ts = Models.DelayDMDc.dt; 
        case 'DMDc'
            this_model.xmean = xref;
        case 'SINDYc'
            this_model.ahat = this_model.Xi(:,1:2);         
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
            COSTFUN = @(u) lotkaObjectiveFCN_models(u,xhat,N,xref,uopt(1),this_model,diag(Q),R,Ru,select_model);
            CONSFUN = @(u) lotkaConstraintFCN_models(u,xhat,N,this_model,select_model);
            uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);
            %   uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options);
            
            
            % Integrate system
            xhat = rk4u(@lotkacontrol_discrete,xhat,uopt(1),Ts/1,1,[],p); %10, 2
            xHistory(:,ct+1) = xhat;
            uHistory(:,ct+1) = uopt(1);
            tHistory(:,ct+1) = ct*Ts+Tcontrol;
            
            if strcmp(select_model,'DelayDMDc')
                if  pest.Nxlim-Ndelay+1+ct<=pest.Nxlim
                    pest.udelay = [zeros(1,Ndelay-ct),uHistory(1:ct)];
                    pest.xdelay = zeros(2,length(pest.udelay));
                    pest.xdelay(:,1:Ndelay) = [x(end-Ndelay+1+ct:end,1:2);xHistory(1:2,1:ct)']';
                else
                    pest.udelay = [uHistory(ct-Ndelay+1:ct)];
                    pest.xdelay = zeros(2,length(pest.udelay));
                    pest.xdelay = [xHistory(1:2,ct-Ndelay+1:ct)];
                end
                if size(pest.xdelay,2)~=Ndelay
                    disp('Error in size.')
                    return
                end
                if length(pest.udelay)~=Ndelay
                    disp('Error in size.')
                    return
                end
            end
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
[tU,xU] = ode45(@(t,x)lotkacontrol(t,x,0,a,b,d,g),[Tcontrol:Ts:Tcontrol+Duration],xA(end,1:2),options);   % true model