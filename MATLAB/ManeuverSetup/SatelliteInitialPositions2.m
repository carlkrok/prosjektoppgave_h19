% T = 0; % 
% i = 24.5; % 
% O = 0;
% w = 0;
% rSat_a = rEarth + 309
% rSat_p = rEarth + 605
r0ECI_target = [6983.13660000000;0;0];
v0ECI_target = [0;6.80006439998294;3.09896788638466];
period_target = 5623.82474323939;

% % No inclination
% r0ECI_chaser = [6587.66508650894;-1056.99562901157;-481.700660123752];
% v0ECI_chaser = [1.32637812496519;6.99546836522539;3.18801860378082];
% period_chaser = [5623.82474323939];
% maneuverStartDelay = 2160;
% orbitType = "prograde";
% maneuverTime = 4560; % Seconds
% 
% % With inclination
% % r0ECI_chaser = [ 6740.78083128987, 520.711196355516, 276.86705427396 ]';
% % v0ECI_chaser = [ -0.66382767817634, 6.80365264836464, 3.61756628288408 ]';
% % period_chaser = 5670.1696304242;
% % maneuverStartDelay = 4440; %300;
% % orbitType = "prograde";
% % maneuverTime = 660; %2460; % Seconds

% T = -10; % 
% i = 22; % 
% O = 0;
% w = 5;
% rSat_a = rEarth + 309
% rSat_p = rEarth + 605
r0ECI_chaser = [6954.22539968667;-564.113386432189;-227.916602463287];
v0ECI_chaser = [0.680137104600880;6.90240127458609;2.78875113617079];
period_chaser = 5623.82474323939;
maneuverStartDelay = 3660;%3000; % minStartDeltaV = 0.0390689366347962
orbitType = "prograde";
maneuverTime = 4800;%5220; % Seconds
%minStartDeltaV =0.0521387832193322