function dy = lotka(t,y,a,b,d,g)
dy = [
    a*y(1) - b*y(1)*y(2);
    d*y(1)*y(2) - g*y(2);
    ];

