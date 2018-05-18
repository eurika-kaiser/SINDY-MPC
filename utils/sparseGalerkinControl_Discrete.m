function ykplus1 = sparseGalerkinControl_Discrete(t,y,u,p)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

polyorder = p.polyorder;
usesine = p.usesine;
ahat = p.ahat;
yPool = poolData([y(p.SelectVars)' u],length(p.SelectVars)+1,polyorder,usesine);
ykplus1 = (yPool*ahat)';