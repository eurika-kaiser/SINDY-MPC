function dy = lotkacontrol(t,y,u,a,b,d,g)
% Copyright 2016, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Sparse Identification of Nonlinear Dynamical Systems with Control"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

dy = [
a*y(1) - b*y(1)*y(2);
d*y(1)*y(2) - g*y(2) + u; %-u*y(2);
];

