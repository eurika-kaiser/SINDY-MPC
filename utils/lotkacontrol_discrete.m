function dy = lotkacontrol_discrete(t,y,u,p)
a = p.a;
b = p.b;
d = p.d;
g = p.g;

dy = [
    a*y(1) - b*y(1)*y(2);
    d*y(1)*y(2) - g*y(2) + u;
    ];

