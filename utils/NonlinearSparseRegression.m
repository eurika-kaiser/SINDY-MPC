function Xi = NonlinearSparseRegression(x,u,dt,options_method,lambda)

Ninputs = min(size(u));
%% Compute Derivative
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
dx(:,size(x,2)+1:size(x,2)+Ninputs) = repmat(0*dx(:,size(x,2)),[1,Ninputs]);

n = size(dx,2);

%% Sparse regression
clear Theta Xi
Theta = poolData(xaug,n,options_method.order,options_method.usesine);
Xi = sparsifyDynamics(Theta,dx,lambda,n);

