
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


% CANDIDATE 2!!!!
% With inclination
% T = -10; 
% i = 100; 
% O = 0;
% w = 5;
% rSat_a = rEarth + 637
% rSat_p = rEarth + 620
r0ECI_chaser = [6971.63491076541;105.914799473061;-600.672676680026];
v0ECI_chaser = [0.656573334996848;-1.30633753434018;7.40860830878767];
period_chaser = [5836.80750425860];
%maneuverStartDelay = 240; % totDeltaV
orbitType = "retrograde";
%maneuverTime = 5460; % totDeltaV - start 0.126
maneuverStartDelay = 420; % startDeltaV 0.0384
maneuverTime = 5400; % startDeltaV

% % CANDIDATE 1
% % With inclination
% % T = 5; 
% % i = 96; 
% % O = 0;
% % w = 0;
% % rSat_a = rEarth + 620
% % rSat_p = rEarth + 637
% r0ECI_chaser = [6988.40958720094;-63.9093938930249;608.057265546706];
% v0ECI_chaser = [-0.657370814907068;-0.784447727793098;7.46352157751033];
% period_chaser = 5836.8075042586;
% %maneuverStartDelay = 60; %1860;  totDV
% orbitType = "retrograde";
% %maneuverTime = 5340;% 3420; % totDV [Seconds] - start 0.1074
% maneuverStartDelay = 720; % startDeltaV 0.0646
% maneuverTime = 5100; % startDeltaV

% % With inclination
% % T = -5; 
% % i = 100; 
% % O = 0;
% % w = -5;
% % rSat_a = rEarth + 637
% % rSat_p = rEarth + 620
% r0ECI_chaser = [6891.85095694499;211.020611828149;-1196.75735942559];
% v0ECI_chaser = [1.31053612089156;-1.29142361236677;7.32402725424072];
% period_chaser = [5836.80750425860];
% maneuverStartDelay = 2340; 
% orbitType = "retrograde";
% maneuverTime = 960; % Seconds

% % CANDIDATE NO
% % With inclination
% % T = 0; 
% % i = 100; 
% % O = 0;
% % w = -5;
% % rSat_a = rEarth + 637
% % rSat_p = rEarth + 620
% r0ECI_chaser = [6971.50657744179;105.912849801512;-600.661619543215];
% v0ECI_chaser = [0.658168294817288;-1.30633753434018;7.40860830878767];
% period_chaser = [5836.80750425860];
% maneuverStartDelay = 180; 
% orbitType = "retrograde";
% maneuverTime = 5280; % Seconds
