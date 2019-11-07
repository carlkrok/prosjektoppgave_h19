
function [Eph, stats] = ephemeris_v3(Y0, N_Step, Step)

options = rdpset('RelTol',1e-13,'AbsTol',1e-16);
[t,yout, stats] = radau(@Accel_v3,(0:Step:N_Step*Step),Y0,options);
%[t,yout] = ode15s(@Accel_v3,(0:Step:N_Step*Step),Y0);
Eph(:,1) = t;
Eph(:,2:7) = yout;

