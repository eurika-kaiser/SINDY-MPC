function dy = lotkacontrol_discrete(t,y,u,p)
% Copyright 2016, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Sparse Identification of Nonlinear Dynamical Systems with Control"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

a = p.a;
b = p.b;
d = p.d;
g = p.g;

dy = [
a*y(1) - b*y(1)*y(2);
d*y(1)*y(2) - g*y(2) + u;%-u.^2;
];

