function x = rk4u(v,x,u,h,n,t,p)

% RK4U   Runge-Kutta scheme of order 4 for control system
%   rk4u(v,X,U,h,n) performs n steps of the scheme for the vector field v
%   using stepsize h on each row of the matrix X
%
%   v(X,U) maps an (m x d)-matrix X and an (m x p)-matrix U
%          to an (m x d)-matrix 
for i = 1:n
    k1 = v(t,x,u,p); 
    k2 = v(t,x + h/2*k1,u,p);
    k3 = v(t,x + h/2*k2,u,p);
    k4 = v(t,x + h*k3,u,p);
    x = x + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end
