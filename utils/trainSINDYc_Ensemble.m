
if eps == 0
    DERIV_NOISE = 0;
end

%% Compute Derivative
Nt = size(x,1);
if any(strcmp(InputSignalType,{'sine2', 'chirp','prbs', 'sphs','mixed', 'noise','unforced'})==1) && DERIV_NOISE == 0%eps==0 % eps~=0
    % compute derivative using fourth order central difference
    % use TVRegDiff if more error
    dx = zeros(Nt-5,Nvar+1,Nic);
    for iIC = 1:Nic
    dx_tmp = zeros(Nt-5,3);
    for i=3:Nt-3
        for k=1:size(x,2)
            dx_tmp(i-2,k) = (1/(12*dt))*(-x(i+2,k,iIC)+8*x(i+1,k,iIC)-8*x(i-1,k,iIC)+x(i-2,k,iIC));
        end
    end
    dx(:,1:5,iIC) = dx_tmp;
    end
    
    % concatenate
    xnew = [x(3:end-3,:,:)];
    unew = [u(3:end-3,:,:)];

    disp(['DONE: computed derivative'])
else %excitation using Gaussian white noise

end

%% Reshape
% uensemble = [u;uv];
% uensemble = u;
X =[];
DX = [];
uensemble = [];
for iIC = 1:Nic
    X = [X; xnew(:,:,iIC)];
    DX = [DX; dx(:,:,iIC)];
    uensemble = [uensemble; unew(:,:,iIC)];
end
% xaug = [X repmat(uensemble(3:end-3,:),[Nic 1])];
xaug = [X uensemble];
DX(:,Nvar+1:size(xaug,2)) = repmat(0*DX(:,Nvar),[1 size(u,2)]);
    
M = size(xaug,1);

n = size(DX,2)-1;
%% Sparse regression
clear Theta Xi
Theta = poolData(xaug,n,polyorder,usesine);
Theta_norm = zeros(size(Theta,2),1);
for i = 1:size(Theta,2)
   Theta_norm(i) = norm(Theta(:,i));
   Theta(:,i) = Theta(:,i)./Theta_norm(i);
end
m = size(Theta,2);

if exist('lambda_vec') == 1
    Xi = sparsifyDynamicsIndependent(Theta,DX,lambda_vec,n);
else
    Xi = sparsifyDynamics(Theta,DX,lambda,n);
end


if n == 3
    str_vars = {'x','y','u'};
elseif n == 4
    str_vars = {'x','y','z','u'};
elseif n == 5
    str_vars = {'x1','x2','x3','x4','x5'};      
elseif n == 6
    str_vars = {'x1','x2','x3','x4','x5','u'};       
end

for i = 1:size(Theta,2)
   Xi(i,:) = Xi(i,:)./Theta_norm(i);
end

yout = poolDataLIST(str_vars,Xi,n,polyorder,usesine);