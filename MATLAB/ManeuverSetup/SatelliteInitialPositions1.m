
r0ECI_target = [ 6993.1366, 0, 0 ]';
v0ECI_target = [ 0, -1.051360468027, 7.48081844151098 ]';
period_target = 5830.56082013569;

% % Near
% % r0ECI_chaser = [ 6992.07280267906, 16.9856712482927, -120.859330914474 ]';
% % v0ECI_chaser = [ 0.131681493347382, -1.05120053489446, 7.47968045814213 ]';
% % period_chaser = 5830.56082013569;
% % maneuverStartDelay = 300;
% % orbitType = "retrograde";
% % maneuverTime = 5460; % Seconds

% % Far
% r0ECI_chaser = [ 6966.55774794289, 84.8252850570378, -603.563264987457 ]';
% v0ECI_chaser = [ 0.657605494122236, -1.04736457505333, 7.4523861856323 ]';
% period_chaser = 5830.56082013569;
% maneuverStartDelay = 2940;%300;
% orbitType = "retrograde"; %"prograde";%
% maneuverTime = 1980;%5400; % Seconds


% % With inclination
% % T = 355; 
% % i = 96; 
% % O = 0;
% % w = 0;
% % rSat_a = rEarth + 620
% % rSat_p = rEarth + 637
% r0ECI_chaser = [6988.40958720094;63.9093938930249;-608.057265546706];
% v0ECI_chaser = [0.657370814907068;-0.784447727793098;7.46352157751033];
% period_chaser = [5836.80750425860];
% maneuverStartDelay = 4620; 
% orbitType = "retrograde";
% maneuverTime = 660; % Seconds

% CANDIDATE 1
% With inclination
% T = 5; 
% i = 96; 
% O = 0;
% w = 0;
% rSat_a = rEarth + 620
% rSat_p = rEarth + 637
r0ECI_chaser = [6988.40958720094;-63.9093938930249;608.057265546706];
v0ECI_chaser = [-0.657370814907068;-0.784447727793098;7.46352157751033];
period_chaser = 5836.8075042586;
maneuverStartDelay = 1860; 
orbitType = "retrograde";
maneuverTime = 3420; % Seconds

% % With inclination
% % T = 5; 
% % i = 96; 
% % O = 0;
% % w = -5;
% % rSat_a = rEarth + 620
% % rSat_p = rEarth + 637
% r0ECI_chaser = [7015.10417651042;5.04017020625450e-15;-9.26196667542843e-14];
% v0ECI_chaser = [-0.000797479910219520;-0.787451483672596;7.49210040567938];
% period_chaser = 5836.8075042586;
% maneuverStartDelay = 2700; 
% orbitType = "retrograde";
% maneuverTime = 120; % Seconds

% % With inclination
% % T = 0; 
% % i = 100; 
% % O = 0;
% % w = -5;
% % rSat_a = rEarth + 620
% % rSat_p = rEarth + 637
% r0ECI_chaser = [6988.44188730935;106.170135211835;-602.120757613231];
% v0ECI_chaser = [0.656573334996848;-1.30317184572283;7.39065479648552];
% period_chaser = [5836.8075042586];
% maneuverStartDelay = 300; 
% orbitType = "retrograde";
% maneuverTime = 4620; % Seconds
