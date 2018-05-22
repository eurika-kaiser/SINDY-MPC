function [sysmodel_DMDc,U,Up] = DelayDMDc_MV(Hx,Hu,r1,r2,dt,numOutputs,numInputs,numVar,X,Xp)
stackmax = size(Hx,1);
stackmax = stackmax/numVar;

if nargin < 10
    X = Hx(:,1:end-1);
    Xp = Hx(:,2:end);
end
Nt = size(X,2);
Gamma = Hu(:,1:Nt);
Omega = [X; Gamma];
[U,S,V] = svd(Omega,'econ'); 
idx = find(diag(S)<1e-10,1,'first'); 
if isempty(idx)==1
    r1 = size(S,1);
else 
    r1 = idx-1;
end
U = U(:,1:r1);S = S(1:r1,1:r1);V = V(:,1:r1);
G = Xp*V*S^(-1)*U';
[Up,Sp,Vp] = svd(Xp,'econ'); 
idx = find(diag(Sp)<1e-10,1,'first'); 
if isempty(idx)==1
    r2 = size(Sp,1);
else 
    r2 = idx-1;
end
Up = Up(:,1:r2);
Ar = G(:,1:stackmax*numVar);
Br = G(:,stackmax*numVar+1:end);
Cr = zeros(numVar,stackmax*numVar); 
if numVar*stackmax == numOutputs
    Cr = eye(numVar*stackmax);
else
    for i = 1:numVar
        Cr(i,i*stackmax) = 1; % 2017-02-21
    end
end

sysmodel_DMDc = ss(Ar,Br,Cr,zeros(numOutputs,numInputs),dt);
