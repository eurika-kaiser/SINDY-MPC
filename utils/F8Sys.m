function dy = F8Sys(t,y,u,p)
% Nonlinear automatic flight control system
% of an F-8 at Mach = 0.85 and an altitude of 30 000 ft.
% Control of angle of attach.
% Input u: tail deflection angle
% States x: angle of attack, pitch angle, pitch rate
% Tracking reference: r(t) = 0.4 (-0.5/(1+exp(t-8)) + 1/(1+exp(t-30)) -0.4)
% Initial state: x(0) = [0.1 0 0], u(0) = 0
% N = 8, Nu = 3, Q=25I, R=0.1I, Fs = 10Hz, -0.5<=du<=0.5, -0.2<=dy<=0.4

try
    dy = [-0.877*y(1) + y(3) - 0.088*y(1)*y(3) + 0.47*y(1)^2 - 0.019*y(2)^2 - y(1)^2*y(3) + 3.846*y(1)^3 - 0.215*u + 0.28*y(1)^2*u + 0.63*u^3 + 0.47*y(1)*u^2;
            y(3);
            -4.208*y(1) - 0.396*y(3) - 0.47*y(1)^2 - 3.564*y(1)^3 - 20.967*u + 46*y(1)*u^2 + 61.4*u^3 + 6.265*y(1)^2*u]; %
catch
    keyboard
end
        

