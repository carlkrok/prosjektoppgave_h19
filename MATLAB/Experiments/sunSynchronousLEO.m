clc;
clear all;

run earthParameters.m;

% Input
Alt = 570:1:580;          % Altitude,Low Earth orbit (LEO)
a   = Alt + rEarth;       % Mean semimajor axis [km]
e   = 0.0;                % Eccentricity
h = a*(1 - e^2);          % [km]
n = (muEarth./a.^3).^0.5; % Mean motion [s-1]
tol = 1e-18;              % Error tolerance
% Initial guess for the orbital inclination
i0 = 180/pi*acos(-2/3*(h/rEarth).^2*wsEarth./(n*J2Earth));
err = 1e1;
while(err >= tol )
    % J2 perturbed mean motion
    np  = n.*(1 + 1.5*J2Earth*(rEarth./h).^2.*(1 - e^2)^0.5.*(1 - 3/2*sind(i0).^2));
    i = 180/pi*acos(-2/3*(h/rEarth).^2*wsEarth./(np*J2Earth));
    err = abs(i - i0);
    i0 = i;
end
plot(Alt,i,'.b');
grid on;hold on;
xlabel('Altitude,Low Earth orbit (LEO)');
ylabel('Mean orbital inclination');
title('Sun-Synchronous Circular Orbit,Inclination vs Altitude(LEO,J2 perturbed)');
hold off;