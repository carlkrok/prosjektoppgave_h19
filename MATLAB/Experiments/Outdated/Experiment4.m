%% Clear Workspace
clc
clear all
format long g
close all


%% Load Parameters

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC deltaVManeuverStart_chaser

run earthParameters; 
run satelliteParameters;

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0,...
                  'Thrust',0, 'stepCounter', 0, 'thrustStartTime', 0, 'thrustEndTime', 0, 'velocityChangeLVLH', 0, 'thrustLVLHAcceleration', 0);

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 2000;

SAT_Const
constants
load DE436Coeff.mat
PC = DE436Coeff; % JPL Planetary Ephemeris DECEMBER 14 1949 - DECEMBER 21 2149

run loadHPOPFileData

% epoch state (Envisat)
fid = fopen('InitialState_OrbitNTNU.txt','r');
tline = fgetl(fid);
year = str2num(tline(1:4));
mon = str2num(tline(6:7));
day = str2num(tline(9:10));
hour = str2num(tline(12:13));
min = str2num(tline(15:16));
sec = str2num(tline(18:23));
tline = fgetl(fid);
Y0(1) = str2num(tline);
tline = fgetl(fid);
Y0(2) = str2num(tline);
tline = fgetl(fid);
Y0(3) = str2num(tline);
tline = fgetl(fid);
Y0(4) = str2num(tline);
tline = fgetl(fid);
Y0(5) = str2num(tline);
tline = fgetl(fid);
Y0(6) = str2num(tline);
tline = fgetl(fid);
AuxParam.area_solar = str2num(tline(49:end));
tline = fgetl(fid);
AuxParam.area_drag = str2num(tline(38:end));
tline = fgetl(fid);
AuxParam.mass = str2num(tline(19:end));
tline = fgetl(fid);
AuxParam.Cr = str2num(tline(5:end));
tline = fgetl(fid);
AuxParam.Cd = str2num(tline(5:end));
fclose(fid);

% epoch
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
%Y0 = ECEF2ECI(Mjd_UTC, Y0); Already in ECI


%% Experiment Setup

maneuverEndTime = 2500;
maneuverStartTime = 300;
Step   = 0.05;   % [s]
N_Step = maneuverEndTime*1/Step; %             
              
Mjd0   = Mjd_UTC;

AuxParam.thrustDuration = 10;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 40;
AuxParam.m       = 40;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.sRad    = 1;
AuxParam.drag    = 1;
AuxParam.SolidEarthTides = 0;
AuxParam.OceanTides = 1;
AuxParam.Relativity = 0;
AuxParam.Thrust = 1;
  


% shorten PC, eopdata, swdata, Cnm, and Snm
num = fix(N_Step*Step/86400)+2;
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

numSamplePointsInitialTrajectorySimple = 200;
numSamplePointsFinalTrajectorySimple = 1000;

%% Monte Carlo Experiment Setup

MCsampleNum = 101;

meanDeviationTimeSetup = 0;
maxDeviationTimeSetup = 5;

MCdeviationTimes = -maxDeviationTimeSetup : (2*maxDeviationTimeSetup)/(MCsampleNum-1) : maxDeviationTimeSetup;

experimentStartTimeIdeal = maneuverStartTime;


%% Initial Orbit Determination

r0ECI_target = [ 6952.13623, 0.0, 0.0 ]';
v0ECI_target = [ 0.0, -0.264590021956698, 7.57686688488419 ]';

r0ECI_chaser = [ Y0(1), Y0(2), Y0(3) ]'./10^3;
v0ECI_chaser = [ Y0(4), Y0(5), Y0(6) ]'./10^3;

orbitPeriod_chaser = orbitPeriod( muEarth, hNorm_chaser, e_chaser );

% A's position at experiment start
%[ rECIExperimentStart_target, vECIExperimentStart_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, experimentStartTimeIdeal, anomalyErrorTolerance, anomalyMaxIterations );
    
% A's position after ideal manouver
[ rECIManouverEnd_target, vECIManouverEnd_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverEndTime, anomalyErrorTolerance, anomalyMaxIterations );
 
% B's ideal position at experiment start
[ rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStartECI_chaser, deltaVEndECI_chaser, vIntersectOrbit_chaser ] = interceptOrbit( rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser, rECIManouverEnd_target, vECIManouverEnd_target, maneuverEndTime - maneuverStartTime, orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

[ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser );
deltaVStartLVLH_chaser = QmatECItoLVLH_chaser * deltaVStartECI_chaser;
deltaVEndLVLH_chaser = QmatECItoLVLH_chaser * deltaVEndECI_chaser;

AuxParam.velocityChangeLVLH = deltaVStartLVLH_chaser*1000;
AuxParam.thrustLVLHAcceleration = AuxParam.velocityChangeLVLH./AuxParam.thrustDuration;


%% Monte Carlo Experimetns



%%%%%%%%%%%%%  SIMPLE MODEL

% B's position at first iteration start
[ rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + MCdeviationTimes( 1 ), anomalyErrorTolerance, anomalyMaxIterations );

[ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser );
QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;


%%%%%%%%%%%%%  HPOP MODEL
% propagation
AuxParam.thrustLVLHAcceleration = AuxParam.velocityChangeLVLH./AuxParam.thrustDuration;
AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + MCdeviationTimes( 1 ))/86400);
[Eph, stats] = ephemeris_Experiment2(Y0, N_Step, Step);

MC_HPOP_PosEnd( 1, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;


for experimentIndex = 2 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );
    
    
    %%%%%%%%%%%%%  SIMPLE MODEL
    % B's position at maneuver start
    [ rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );
 
    [ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser );
    QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
    deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;
        
    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    AuxParam.thrustLVLHAcceleration = AuxParam.velocityChangeLVLH./AuxParam.thrustDuration;
    AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + thisDeviationTime)/86400);
    [thisEph, stats] = ephemeris_Experiment2(Y0, N_Step, Step);
    MC_HPOP_PosEnd( experimentIndex, : ) = [thisEph(N_Step+1, 2), thisEph(N_Step+1, 3), thisEph(N_Step+1, 4)]./10^3;

   
end

%%

figure(1)
hold on
grid on
title('End Position Error HPOP')
for plotIndex = 1 : MCsampleNum
    plot(MCdeviationTimes( plotIndex ), norm(MC_HPOP_PosEnd( plotIndex, : )' - rECIManouverEnd_target), '*')
end
xlabel('Delay Time [s]')
ylabel('Distance [km]')
hold off

figure(2)
hold on
grid on
title('End Position HPOP')
plot3(0,0,0,'k.')
for plotIndex = 1 : MCsampleNum
	plot3((MC_HPOP_PosEnd( plotIndex, 1 ) - rECIManouverEnd_target(1)), (MC_HPOP_PosEnd( plotIndex, 2 ) - rECIManouverEnd_target(2)), (MC_HPOP_PosEnd( plotIndex, 3 ) - rECIManouverEnd_target(3)), '*')
    if mod((plotIndex-1), 10) == 0
        plot3([0 (MC_HPOP_PosEnd( plotIndex, 1 ) - rECIManouverEnd_target(1))],[0 (MC_HPOP_PosEnd( plotIndex, 2 ) - rECIManouverEnd_target(2))],[0 (MC_HPOP_PosEnd( plotIndex, 3 ) - rECIManouverEnd_target(3))])
    end
end
hold off

