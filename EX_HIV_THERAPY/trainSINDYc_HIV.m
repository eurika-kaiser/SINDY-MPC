
% Compute  derivative using fourth order central difference
% use TVRegDiff if more error
dx = zeros(length(x)-5,NSindy);
for i=3:length(x)-3
    for k=1:NSindy
        dx(i-2,k) = (1/(12*dt))*(-x(i+2,SelectVars(k))+8*x(i+1,SelectVars(k))-8*x(i-1,SelectVars(k))+x(i-2,SelectVars(k)));
    end
end
% concatenate
xaug = [x(3:end-3,SelectVars) u(3:end-3,:)];
dx(:,NSindy+1) = 0*dx(:,1);


n = size(dx,2);

% Sparse regression
clear Theta Xi
Theta = poolData(xaug,n,polyorder,usesine);
Theta_norm = zeros(size(Theta,2),1); %zeros(size(Theta,2),1);ones
for i = 1:size(Theta,2)
    Theta_norm(i) = norm(Theta(:,i));
    Theta(:,i) = Theta(:,i)./Theta_norm(i);
end
m = size(Theta,2);

if exist('lambda_vec') == 1
    Xi = sparsifyDynamicsIndependent(Theta,dx,lambda_vec,n-1);
else
    Xi = sparsifyDynamics(Theta,dx,lambda,n-1);
end

% count_lambda = 0;
% while all(Xi(:)==0)==1 && count_lambda < 5
%     lambda = lambda/10
%     Xi = sparsifyDynamics(Theta,dx,lambda,n);
%     count_lambda = count_lambda + 1;
% end

str_vars = {'x1','x2','x3','x4','x5','u'};
str_vars = str_vars([SelectVars,end]);

for i = 1:size(Theta,2)
    Xi(i,:) = Xi(i,:)./Theta_norm(i);
end

yout = poolDataLIST(str_vars,Xi,n,polyorder,usesine);