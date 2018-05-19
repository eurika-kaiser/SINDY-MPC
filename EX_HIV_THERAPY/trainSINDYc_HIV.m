%% Compute time derivative
% Compute  derivative using fourth order central difference
% use TVRegDiff if more error
dx = zeros(length(x)-5,NSindy);
for i=3:length(x)-3
    for k=1:NSindy
        dx(i-2,k) = (1/(12*dt))*(-x(i+2,SelectVars(k))+8*x(i+1,SelectVars(k))-8*x(i-1,SelectVars(k))+x(i-2,SelectVars(k)));
    end
end

% Concatenate states, input, and time derivatives
xaug = [x(3:end-3,SelectVars) u(3:end-3,:)];
dx(:,NSindy+1) = 0*dx(:,1);
n = size(dx,2);

%% Sparse regression
% Construct library
clear Theta Xi
Theta = poolData(xaug,n,polyorder,usesine);

% Normalize library columns
Theta_norm = zeros(size(Theta,2),1); 
for i = 1:size(Theta,2)
    Theta_norm(i) = norm(Theta(:,i));
    Theta(:,i) = Theta(:,i)./Theta_norm(i);
end
m = size(Theta,2);

% Estimate parameters
if exist('lambda_vec') == 1 % choose different penalization for each variable
    Xi = sparsifyDynamicsIndependent(Theta,dx,lambda_vec,n-1);
else % same penalization for all variables
    Xi = sparsifyDynamics(Theta,dx,lambda,n-1);
end

% Rescale identified coefficients to account for normalization
for i = 1:size(Theta,2)
    Xi(i,:) = Xi(i,:)./Theta_norm(i);
end

% Print coefficients
str_vars = {'x1','x2','x3','x4','x5','u'};
str_vars = str_vars([SelectVars,end]);
yout = poolDataLIST(str_vars,Xi,n,polyorder,usesine);