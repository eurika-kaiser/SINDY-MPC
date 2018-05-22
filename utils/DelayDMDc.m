function [sysmodel_DMDc,U,Up] = DelayDMDc(Hx,Hu,r1,r2,dt,numOutputs,numInputs)
stackmax = size(Hx,1);

X = Hx(:,1:end-1);
Xp = Hx(:,2:end);
Nt = size(X,2);
Gamma = Hu(:,1:Nt);
Omega = [X; Gamma];
[U,S,V] = svd(Omega,'econ'); 
U = U(:,1:r1);S = S(1:r1,1:r1);V = V(:,1:r1);
G = Xp*V*S^(-1)*U';
[Up,Sp,Vp] = svd(Xp,'econ'); Up = Up(:,1:r2);
Ar = Up'*G(:,1:stackmax)*Up;
Br = Up'*G(:,stackmax+1:end);
% Cr = Up(1,:);
%Cr = Up(end,:);
Cr = zeros(1,stackmax); Cr(stackmax) = 1; % 2017-02-21
Cr = Cr*Up; % == Up(end,:)

% Ar = U(1:stackmax,1:r2)'*G(:,1:stackmax)*U(1:stackmax,1:r2);
% Br = U(1:stackmax,1:r2)'*G(:,stackmax+1:end);
% Cr = U(stackmax,1:r2);
sysmodel_DMDc = ss(Ar,Br,Cr,zeros(numOutputs,numInputs),dt);
