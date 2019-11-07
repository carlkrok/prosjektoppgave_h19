%% Clear Workspace
clc
clear all
format long g
close all


%% Load Parameters

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC 

run earthParametersHPOP; 

r0ECI_target = [ 6990.1366, 0.0, 0.0 ]';
v0ECI_target = [ 0.0, -1.04884571201957, 7.50091780930118 ]';

r0ECI_chaser = [ 6952.1366, 0.0, 0.0 ]';
v0ECI_chaser = [ 0.0, -0.976461903728063, 7.51834019814233 ]';

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0,...
                  'Thrust',0, 'stepCounter', 0, 'velocityChangeLVLH', 0, 'thrustLVLHAcceleration', 0);

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 2000;

SAT_Const
constants
load DE436Coeff.mat
PC = DE436Coeff; % JPL Planetary Ephemeris DECEMBER 14 1949 - DECEMBER 21 2149


% read Earth gravity field coefficients
Cnm = zeros(181,181); % Geopotential coefficient (cos)
Snm = zeros(181,181); % Geopotential coefficient (sin)
fid = fopen('GGM03S.txt','r'); % GRACE Gravity Model 03 (GGM03), GGM03S - complete to harmonic degree 180
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

% read Earth orientation parameters 1999 - 2018
fid = fopen('eop19990101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

% read space weather data 1999 - 2017
fid = fopen('sw19990101.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
swdata = fscanf(fid,'%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4f %2i %4i %6f %2i %6f %6f %6f %6f %6f',[33 inf]);
fclose(fid);

% read space weather data 1997 - 2016 Solar data
fid = fopen('SOLFSMY.txt','r'); % SOLFSMY_2019
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% READ GEOMAGNETIC STORM DTC VALUE
fid = fopen('DTCFILE.txt','r'); % DTCFILE_2019
%  ------------------------------------------------------------------------
% | DTC YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);

% read space weather data
fid = fopen('SOLRESAP.txt','r'); % SOLRESAP_2019
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[12 inf]);
fclose(fid);


% epoch state
year = 2002;
mon = 04;
day = 24;
hour = 00;
min = 00;
sec = 00;
Y0 = [ r0ECI_chaser', v0ECI_chaser'].*10^3;
AuxParam.area_solar = 0.2;
AuxParam.area_drag = 0.1;
AuxParam.mass = 2.0;
AuxParam.Cr = 1.0;
AuxParam.Cd = 2.2;

% epoch
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);

%% Experiment Setup

maneuverEndTime = 2500;
maneuverStartTime = 300;
Step   = 0.5;   % [s]
N_Step = round(maneuverEndTime*1/Step);
N_Step_Initial = round(maneuverStartTime*1/Step); % 
              
Mjd0   = Mjd_UTC;

AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 0;
AuxParam.m       = 0;
AuxParam.sun     = 0;
AuxParam.moon    = 0;
AuxParam.planets = 0;
AuxParam.sRad    = 0;
AuxParam.drag    = 0;
AuxParam.SolidEarthTides = 0;
AuxParam.OceanTides = 0;
AuxParam.Relativity = 0;
AuxParam.Thrust = 0;

% shorten PC, eopdata, swdata, Cnm, and Snm
num = fix(N_Step*Step/const.DAYSEC)+2;
JD = Mjd_UTC+2400000.5;
i = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
PC = PC(i:i+num,:);
mjd = (floor(Mjd_UTC));
i = find(mjd==eopdata(4,:),1,'first');
eopdata = eopdata(:,i:i+num);
i = find((year==swdata(1,:)) & (mon==swdata(2,:)) & (day==swdata(3,:)),1,'first');
swdata = swdata(:,i-3:i+num);
Cnm = Cnm(1:AuxParam.n+1,1:AuxParam.n+1);
Snm = Snm(1:AuxParam.n+1,1:AuxParam.n+1);

orbitType_chaser = "retrograde";


%% Initial Orbit Determination

% A's position after ideal manouver
%[ rECIManouverEnd_target, vECIManouverEnd_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverEndTime, anomalyErrorTolerance, anomalyMaxIterations );
targetY0 = [ r0ECI_target', v0ECI_target' ].*10^3;
[tagetEph] = ephemeris_v3(targetY0, N_Step, Step);
rECIManouverEnd_target = [ tagetEph( N_Step+1, 2 ); tagetEph( N_Step+1, 3 ); tagetEph( N_Step+1, 4 ) ]./10^3;
vECIManouverEnd_target = [ tagetEph( N_Step+1, 5 ); tagetEph( N_Step+1, 6 ); tagetEph( N_Step+1, 7 ) ]./10^3;

% B's ideal position at experiment start
%[ rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, maneuverStartTime, anomalyErrorTolerance, anomalyMaxIterations );
chaserY0 = [ r0ECI_chaser', v0ECI_chaser' ].*10^3;
[chaserEph] = ephemeris_v3(chaserY0, N_Step_Initial, Step);
rECIManeuverStartIdeal_chaser = [ chaserEph( N_Step_Initial+1, 2 ); chaserEph( N_Step_Initial+1, 3 ); chaserEph( N_Step_Initial+1, 4 ) ]./10^3;
vECIManeuverStartIdeal_chaser = [ chaserEph( N_Step_Initial+1, 5 ); chaserEph( N_Step_Initial+1, 6 ); chaserEph( N_Step_Initial+1, 7 ) ]./10^3;

% Required velocity change satellite B
[ deltaVStartECI_chaser, deltaVEndECI_chaser, vIntersectOrbit_chaser ] = interceptOrbit( rECIManeuverStartIdeal_chaser, vECIManeuverStartIdeal_chaser, rECIManouverEnd_target, vECIManouverEnd_target, maneuverEndTime - maneuverStartTime, orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

[ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartIdeal_chaser, vECIManeuverStartIdeal_chaser );
deltaVStartLVLH_chaser = QmatECItoLVLH_chaser * deltaVStartECI_chaser;

AuxParam.velocityChangeLVLH = deltaVStartLVLH_chaser*1000;


%% Experimetns

thrustDurationVec = 1 : 1 : 15;
numThrustDuration = length(thrustDurationVec);
HPOP_PosEnd = zeros(numThrustDuration,3);


%%%%%%%%%%%%%  SIMPLE MODEL


% B's position at first iteration end
[ rECIExperimentEndSimple_chaser, vECIExperimentEndSimple_chaser ] = nextStateTimeStep( muEarth, rECIManeuverStartIdeal_chaser, vECIManeuverStartIdeal_chaser + deltaVStartECI_chaser, maneuverEndTime - maneuverStartTime, anomalyErrorTolerance, anomalyMaxIterations );


for experimentIndex = 1 : numThrustDuration
    
    %%%%%%%%%%%%%  HPOP MODEL
    
    AuxParam.thrustDuration = thrustDurationVec( experimentIndex );
    AuxParam.thrustLVLHAcceleration = AuxParam.velocityChangeLVLH./AuxParam.thrustDuration;

    N_Step_Thrust = round(AuxParam.thrustDuration*1/Step); %  
    N_Step_Final = round(maneuverEndTime*1/Step) - N_Step_Initial - N_Step_Thrust; %  
    
    % propagation
    
    AuxParam.Thrust = 0;
    [initialEph] = ephemeris_v3(Y0, N_Step_Initial, Step);
    currY = [ initialEph(N_Step_Initial+1, 2), initialEph(N_Step_Initial+1, 3), initialEph(N_Step_Initial+1, 4), initialEph(N_Step_Initial+1, 5), initialEph(N_Step_Initial+1, 6), initialEph(N_Step_Initial+1, 7) ];
    
    AuxParam.Thrust = 1;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime)/const.DAYSEC);
    [thrustEph] = ephemeris_v3(currY, N_Step_Thrust, Step);
    currY = [ thrustEph(N_Step_Thrust+1, 2), thrustEph(N_Step_Thrust+1, 3), thrustEph(N_Step_Thrust+1, 4), thrustEph(N_Step_Thrust+1, 5), thrustEph(N_Step_Thrust+1, 6), thrustEph(N_Step_Thrust+1, 7) ];
    
    AuxParam.Thrust = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime)/const.DAYSEC) + ((AuxParam.thrustDuration)/const.DAYSEC);
    [finalEph] = ephemeris_v3(currY, N_Step_Final, Step);

    HPOP_PosEnd( experimentIndex, : ) = [finalEph(N_Step_Final+1, 2), finalEph(N_Step_Final+1, 3), finalEph(N_Step_Final+1, 4)]./10^3;

   
end

%%

[ QmatECItoLVLH_target ] = ECIToLVLH( rECIManouverEnd_target, vECIManouverEnd_target );

endPosLVLH = zeros(numThrustDuration,3)
for posIter = 1:numThrustDuration
    thisPos = HPOP_PosEnd( posIter, : )'- rECIManouverEnd_target;
    thisPos = QmatECItoLVLH_target * thisPos;
    endPosLVLH( posIter, : ) = thisPos';
end

%%

figure(1)
hold on
grid on
title('End Position Error HPOP')
for plotIndex = 1 : numThrustDuration
    plot(thrustDurationVec( plotIndex ), norm(HPOP_PosEnd( plotIndex, : )' - rECIManouverEnd_target), '*')
end
xlabel('Maneuver Time [s]')
ylabel('Distance [km]')
hold off


figure(2)
hold on
grid on
title('End Position ECI HPOP')
plot3(0,0,0,'k.')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
for plotIndex = 1 : numThrustDuration
	plot3((HPOP_PosEnd( plotIndex, 1 ) - rECIManouverEnd_target(1)), (HPOP_PosEnd( plotIndex, 2 ) - rECIManouverEnd_target(2)), (HPOP_PosEnd( plotIndex, 3 ) - rECIManouverEnd_target(3)), '*')
end
hold off


figure(3)
hold on
grid on
title('End Position LVLH HPOP')
plot3(0,0,0,'k.')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
for plotIndex = 1 : numThrustDuration
	plot3(endPosLVLH( plotIndex, 1 ), endPosLVLH( plotIndex, 2 ), endPosLVLH( plotIndex, 3 ), '*')
end
hold off
