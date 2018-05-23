
if eps == 0
    DERIV_NOISE = 0;
end


if size(x,1)==size(u,2)
    u = u';
end

% Compute Derivative
if any(strcmp(InputSignalType,{'sine2', 'sine3', 'chirp','prbs', 'sphs','mixed', 'noise','unforced'})==1) && DERIV_NOISE == 0%eps==0 % eps~=0
    % compute derivative using fourth order central difference
    % use TVRegDiff if more error
    dx = zeros(length(x)-5,3);
    for i=3:length(x)-3
        for k=1:size(x,2)
            dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k));
        end
    end
    % concatenate
    xaug = [x(3:end-3,:) u(3:end-3,:)];
    dx(:,size(x,2)+1) = 0*dx(:,size(x,2));
else %excitation using Gaussian white noise
    %xclean = x;
    dxclean = zeros(size(xclean));
    % Compute clean derivative  (just for comparison!)
    for i=1:length(x)
        dxclean(i,:) = LorenzSys(0,xclean(i,:),u(i),p);
    end
    xclean = [xclean u];
    dxclean(:,4) = 0*dxclean(:,3);

    %  Total Variation Regularized Differentiation
    clear dx xt
    for i = 1:size(x,2)
        dx(:,i) = TVRegDiff( x(:,i), 20, .00002, [], 'small', 1e12, dt, 1, 1 ); %.00002
    end
    dx = dx(2:end,:);
    
    for i = 1:size(x,2)
        xt(:,i) = cumsum(dx(:,i))*dt;
        xt(:,i) = xt(:,i) - (mean(xt(50:end-50,i)) - mean(x(50:end-50,i)));
    end
    xt = xt(50:end-51,:);
    dx = dx(50:end-51,:);  % trim off ends (overly conservative)
    xaug = [xt u(50:end-51)];
    dx(:,size(x,2)+1) = 0*dx(:,size(x,2));
    close all
end

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


if n == 3
    str_vars = {'x','y','u'};
elseif n == 4
    str_vars = {'x','y','z','u'};
elseif n == 6
    str_vars = {'x1','x2','x3','x4','x5','u'};    
end

for i = 1:size(Theta,2)
   Xi(i,:) = Xi(i,:)./Theta_norm(i);
end

yout = poolDataLIST(str_vars,Xi,n,polyorder,usesine);