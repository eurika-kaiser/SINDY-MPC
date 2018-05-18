function u = forcing_sine2_positive(x,t)

u = [(0.4*(sin(2*pi*0.2*t).*sin(2*pi*0.05*t)))+0.3];
u(u<0) = 0;
% if u<0
%     u = 0;
% end

