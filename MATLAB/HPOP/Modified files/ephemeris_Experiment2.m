
function [Eph, stats] = ephemeris_Experiment2(Y0, N_Step, Step)

options = rdpset('RelTol',1e-13,'AbsTol',1e-16);
[t,yout, stats] = radau(@Accel_Experiment2,(0:Step:N_Step*Step),Y0,options);
Eph(:,1) = t;
Eph(:,2:7) = yout;

