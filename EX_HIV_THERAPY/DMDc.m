function [sysmodel_DMDc,U,Up] = DMDc(X1,X2,U,dt)
numVar = size(X1,1); numOutputs = numVar; numInputs = size(U,1);
Gamma = U;
Omega = [X1; Gamma];
[U,S,V] = svd(Omega,'econ'); 
G = X2*V*S^(-1)*U';
[Up,Sp,Vp] = svd(X2,'econ'); 
Ar = G(:,1:numVar);
Br = G(:,numVar+1:end);
Cr = eye(numVar); 
sysmodel_DMDc = ss(Ar,Br,Cr,zeros(numOutputs,numInputs),dt);
