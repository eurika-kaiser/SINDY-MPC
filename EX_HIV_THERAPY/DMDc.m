function [sysmodel_DMDc,U,Up] = DMDc(X1,X2,U,dt)
numOutputs = size(X1,1); numInputs = size(U,1);
Omega = [X1; U];
[U,S,V] = svd(Omega,'econ'); 
G = X2*V*S^(-1)*U';

% Needed if reduced-order projection is seeked
[Up,Sp,Vp] = svd(X2,'econ'); 

% Construct system matrices
Ar = G(:,1:numOutputs);
Br = G(:,numOutputs+1:end);
Cr = eye(numOutputs); 

% Discrete-time state- space model
sysmodel_DMDc = ss(Ar,Br,Cr,zeros(numOutputs,numInputs),dt);
