function [t,y,ie] = integrateODE(tspan, x0, solver, time_limit,forcing)
% Wrapper to terminate integration when it doesn't converge (when it exceeds duration time_limit)

n = length(x0);
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n),'Events', @(t,y)timeoutFun(t,y,time_limit));
start_time = tic;
switch solver
    case 'ode45'
        [t, y, te, ye, ie] = ode45(@(t,x) HIVsys_ZURAKOWSKI(t,x, forcing(x,t)), tspan, x0, options);
    case 'ode15s'
        [t, y, te, ye, ie] = ode15s(@(t,x) HIVsys_ZURAKOWSKI(t,x, forcing(x,t)), tspan, x0, options);
    case 'ode23tb'
        [t, y, te, ye, ie] = ode23tb(@(t,x) HIVsys_ZURAKOWSKI(t,x, forcing(x,t)), tspan, x0, options);
end
if isempty(ie)==1
    ie = 0;
end

function [value, isterminal, direction] = timeoutFun(t, y,time_limit)

% time_limit = 100; % set a time limit (in seconds)
value1 = time_limit - toc(start_time);
% value2 = min(y);

% value = min([value1, value2]); 
value = value1;

isterminal = 1;
direction = 0;

end
end



