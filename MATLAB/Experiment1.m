%% Clear Workspace
clc
clear all
format long g
close all


%% Load Parameters

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC

run earthParameters; 
run satelliteParameters;

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0,...
                  'Thrust',0, 'stepCounter', 0, 'thrustStartTime', 0, 'thrustEndTime', 0, 'VelocityChange', 0, 'thrustAcceleration', 0);

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

experimentTime = 2500;
maneuverTime = 300;
AuxParam.thrustDuration = 10;
Step   = 0.1;   % [s]
N_Step = experimentTime*1/Step; % 

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
AuxParam.Thrust = 1;
AuxParam.VelocityChange = 0;              
              
Mjd0   = Mjd_UTC;

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

MCsampleNum = 10;

meanDeviationTimeSetup = 0;
maxDeviationTimeSetup = 2;

MCdeviationTimes = meanDeviationTimeSetup - maxDeviationTimeSetup + ( 2 * maxDeviationTimeSetup * rand( MCsampleNum, 1 ) );

experimentStartTimeIdeal = maneuverTime;


%% Initial Orbit Determination

% r0PQW_target = positionVectorPQW( muEarth, hNorm_target, e_target, T_target );
% v0PQW_target = velocityVectorPQW( muEarth, hNorm_target, e_target, T_target );
% QmatPQWtoECI_target = transformPQWtoECI( i_target, O_target, w_target );
% r0ECI_target = QmatPQWtoECI_target * r0PQW_target;
% v0ECI_target = QmatPQWtoECI_target * v0PQW_target;
r0ECI_target = [ 6952.13623, 0.0, 0.0 ]';
v0ECI_target = [ 0.0, -0.264590021956698, 7.57686688488419 ]';

% r0PQW_chaser = positionVectorPQW( muEarth, hNorm_chaser, e_chaser, T_chaser );
% v0PQW_chaser = velocityVectorPQW( muEarth, hNorm_chaser, e_chaser, T_chaser );
% QmatPQWtoECI_chaser = transformPQWtoECI( i_chaser, O_chaser, w_chaser );
% r0ECI_chaser = QmatPQWtoECI_chaser * r0PQW_chaser;
% v0ECI_chaser = QmatPQWtoECI_chaser * v0PQW_chaser;
r0ECI_chaser = [ Y0(1), Y0(2), Y0(3) ]'./10^3;
v0ECI_chaser = [ Y0(4), Y0(5), Y0(6) ]'./10^3;

orbitPeriod_chaser = orbitPeriod( muEarth, hNorm_chaser, e_chaser );

% A's position at experiment start
%[ rECIExperimentStart_target, vECIExperimentStart_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, experimentStartTimeIdeal, anomalyErrorTolerance, anomalyMaxIterations );
    
% A's position after ideal manouver
[ rECIManouverEnd_target, vECIManouverEnd_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, experimentTime, anomalyErrorTolerance, anomalyMaxIterations );
 
% B's ideal position at experiment start
[ rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStartECI_chaser, deltaVEndECI_chaser, vIntersectOrbit_chaser ] = interceptOrbit( rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser, rECIManouverEnd_target, vECIManouverEnd_target, experimentTime - maneuverTime, orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

[ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser );
deltaVStartLVLH_chaser = QmatECItoLVLH_chaser * deltaVStartECI_chaser;
deltaVEndLVLH_chaser = QmatECItoLVLH_chaser * deltaVEndECI_chaser;


%% Monte Carlo Experimetns

MCSimplePosEnd = zeros( MCsampleNum, 3);
MCSimplePosStart = zeros( MCsampleNum, 3);
MCSimplePosEndMean = zeros(MCsampleNum, 3);
MCSimplePosStartMean = zeros(MCsampleNum, 3);

MCSimpleVelEnd = zeros(MCsampleNum, 3);
MCSimpleVelStart = zeros(MCsampleNum, 3);
MCSimpleVelEndMean = zeros(MCsampleNum, 3);
MCSimpleVelStartMean = zeros(MCsampleNum, 3);

MCSimpleECI_X_Trajectories = zeros(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);
MCSimpleECI_Y_Trajectories = zeros(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);
MCSimpleECI_Z_Trajectories = zeros(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);

MC_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_HPOP_PosEndMean = zeros(MCsampleNum, 3);
MC_HPOP_PosStartMean = zeros(MCsampleNum, 3);

MC_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_HPOP_VelEndMean = zeros(MCsampleNum, 3);
MC_HPOP_VelStartMean = zeros(MCsampleNum, 3);

MC_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);

%%%%%%%%%%%%%  SIMPLE MODEL
% (Debug) B's initial trajectory
[rXECITrajectoryInitial_chaser, rYECITrajectoryInitial_chaser, rZECITrajectoryInitial_chaser, vXECITrajectoryInitial_chaser, vYECITrajectoryInitial_chaser, vZECITrajectoryInitial_chaser, sampleTECITrajectoryInitial_chaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, 1, numSamplePointsInitialTrajectorySimple, muEarth );

% B's position at first iteration start
[ rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + MCdeviationTimes( 1 ), anomalyErrorTolerance, anomalyMaxIterations );

[rXECITrajectoryExperimentStart_chaser, rYECITrajectoryExperimentStart_chaser, rZECITrajectoryExperimentStart_chaser, vXECITrajectoryExperimentStart_chaser, vYECITrajectoryExperimentStart_chaser, vZECITrajectoryExperimentStart_chaser, sampleTECITrajectoryExperimentStart_chaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime + MCdeviationTimes( 1 ), 1, numSamplePointsInitialTrajectorySimple, muEarth );
MCSimpleECI_X_Trajectories(1, 1:numSamplePointsInitialTrajectorySimple) = rXECITrajectoryExperimentStart_chaser;
MCSimpleECI_Y_Trajectories(1, 1:numSamplePointsInitialTrajectorySimple) = rYECITrajectoryExperimentStart_chaser;
MCSimpleECI_Z_Trajectories(1, 1:numSamplePointsInitialTrajectorySimple) = rZECITrajectoryExperimentStart_chaser;

MCSimplePosStart( 1, : ) = rECIExperimentStartSimple_chaser';
MCSimpleVelStart( 1, : ) = vECIExperimentStartSimple_chaser';

MCSimplePosStartMean( 1, : ) = MCSimplePosStart( 1, : );

[ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser );
QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
deltaVExperimentStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;

% B's position at first iteration end
[ rECIExperimentEndSimple_chaser, vECIExperimentEndSimple_chaser ] = nextStateTimeStep( muEarth, rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser + deltaVExperimentStart_chaser, experimentTime - maneuverTime - MCdeviationTimes( 1 ), anomalyErrorTolerance, anomalyMaxIterations );

% (Debug) B's new trajectory
[rXECITrajectoryNew_chaser, rYECITrajectoryNew_chaser, rZECITrajectoryNew_chaser, vXECITrajectoryNew_chaser, vYECITrajectoryNew_chaser, vZECITrajectoryNew_chaser, sampleTECITrajectoryNew_chaser] = ECITrajectory( rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser + deltaVExperimentStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, experimentTime - maneuverTime, 1, numSamplePointsFinalTrajectorySimple, muEarth );

[rXECITrajectoryExperimentEnd_chaser, rYECITrajectoryExperimentEnd_chaser, rZECITrajectoryExperimentEnd_chaser, vXECITrajectoryExperimentEnd_chaser, vYECITrajectoryExperimentEnd_chaser, vZECITrajectoryExperimentEnd_chaser, sampleTECITrajectoryExperimentEnd_chaser] = ECITrajectory( rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser + deltaVExperimentStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, experimentTime - maneuverTime - MCdeviationTimes( 1 ), 1, numSamplePointsFinalTrajectorySimple, muEarth );
MCSimpleECI_X_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rXECITrajectoryExperimentEnd_chaser;
MCSimpleECI_Y_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rYECITrajectoryExperimentEnd_chaser;
MCSimpleECI_Z_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rZECITrajectoryExperimentEnd_chaser;

MCSimplePosEnd( 1, : ) = rECIExperimentEndSimple_chaser';
MCSimpleVelEnd( 1, : ) = vECIExperimentEndSimple_chaser';

MCSimplePosEndMean( 1, : ) = MCSimplePosEnd( 1, : );
MCSimpleVelEndMean( 1, : ) = MCSimpleVelEnd( 1, : );


%%%%%%%%%%%%%  HPOP MODEL
% propagation
AuxParam.thrustAcceleration = (deltaVExperimentStart_chaser*1000)./AuxParam.thrustDuration;
AuxParam.thrustStartTime = Mjd0 + ((maneuverTime + MCdeviationTimes( 1 ))/86400);
[Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
MC_HPOP_PosEnd( 1, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
MC_HPOP_VelEnd( 1, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

MC_HPOP_PosEndMean( 1, : ) = MC_HPOP_PosEnd( 1, : );
MC_HPOP_VelEndMean( 1, : ) = MC_HPOP_VelEnd( 1, : );

MC_HPOP_ECI_X_Trajectories(1, :) = Eph(:, 2)'./10^3;
MC_HPOP_ECI_Y_Trajectories(1, :) = Eph(:, 3)'./10^3;
MC_HPOP_ECI_Z_Trajectories(1, :) = Eph(:, 4)'./10^3;


for experimentIndex = 2 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );
    
    
    %%%%%%%%%%%%%  SIMPLE MODEL
    % B's position at maneuver start
    [ rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    MCSimplePosStart( experimentIndex, : ) = rECIManeuverStartSimple_chaser';
    MCSimpleVelStart( experimentIndex, : ) = vECIManeuverStartSimple_chaser';
    
    [rXECITrajectoryManeuverStart_chaser, rYECITrajectoryManeuverStart_chaser, rZECITrajectoryManeuverStart_chaser, vXECITrajectoryManeuverStart_chaser, vYECITrajectoryManeuverStart_chaser, vZECITrajectoryManeuverStart_chaser, sampleTECITrajectoryManeuverStart_chaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime + MCdeviationTimes( experimentIndex ), 1, numSamplePointsInitialTrajectorySimple, muEarth );
    MCSimpleECI_X_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rXECITrajectoryManeuverStart_chaser;
    MCSimpleECI_Y_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rYECITrajectoryManeuverStart_chaser;
    MCSimpleECI_Z_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rZECITrajectoryManeuverStart_chaser;

    MCSimplePosStartMean( experimentIndex, : ) = mean( MCSimplePosStart( 1 : experimentIndex, : ) );
    
    [ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser );
    QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
    deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;

    % B's position at experiment end
    [ rECIExperimentEndSimple_chaser, vECIExperimentEndSimple_chaser ] = nextStateTimeStep( muEarth, rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser + deltaVManeuverStart_chaser, experimentTime - maneuverTime - thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    MCSimplePosEnd( experimentIndex, : ) = rECIExperimentEndSimple_chaser';
    MCSimpleVelEnd( experimentIndex, : ) = vECIExperimentEndSimple_chaser';
    
    MCSimplePosEndMean( experimentIndex, : ) = mean( MCSimplePosEnd( 1 : experimentIndex, : ) );
    MCSimpleVelEndMean( experimentIndex, : ) = mean( MCSimpleVelEnd( 1 : experimentIndex, : ) );
    
    [rXECITrajectoryManeuverEnd_chaser, rYECITrajectoryManeuverEnd_chaser, rZECITrajectoryManeuverEnd_chaser, vXECITrajectoryManeuverEnd_chaser, vYECITrajectoryManeuverEnd_chaser, vZECITrajectoryManeuverEnd_chaser, sampleTECITrajectoryManeuverEnd_chaser] = ECITrajectory( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser + deltaVManeuverStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, experimentTime - maneuverTime - MCdeviationTimes( experimentIndex ), 1, numSamplePointsFinalTrajectorySimple, muEarth );
    MCSimpleECI_X_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rXECITrajectoryManeuverEnd_chaser;
    MCSimpleECI_Y_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rYECITrajectoryManeuverEnd_chaser;
    MCSimpleECI_Z_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rZECITrajectoryManeuverEnd_chaser;

    
    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    AuxParam.thrustAcceleration = (deltaVExperimentStart_chaser*1000)./AuxParam.thrustDuration;
    AuxParam.thrustStartTime = Mjd0 + ((maneuverTime + thisDeviationTime)/86400);
    [Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
    MC_HPOP_PosEnd( experimentIndex, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
    MC_HPOP_VelEnd( experimentIndex, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

    MC_HPOP_PosEndMean( experimentIndex, : ) = mean(MC_HPOP_PosEnd( 1 : experimentIndex, : ));
    MC_HPOP_VelEndMean( experimentIndex, : ) = mean(MC_HPOP_VelEnd( 1 : experimentIndex, : ));
    
    MC_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph(:, 2)'./10^3;
    MC_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph(:, 3)'./10^3;
    MC_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph(:, 4)'./10^3;
    
end


%% Statistical Analysis

absMeanDeviationEndPosSimple = zeros( MCsampleNum, 1 );
absDeviationEndPosSimple = zeros( MCsampleNum, 1 );
relEndPosSimple_chaser = MCSimplePosEnd;

absMeanDeviationEndPosHPOP = zeros( MCsampleNum, 1 );
absDeviationEndPosHPOP = zeros( MCsampleNum, 1 );
relEndPosHPOP_chaser = MC_HPOP_PosEnd;

for dataIndex = 1 : MCsampleNum
    
    thisMeanEndPosSimple = MCSimplePosEndMean( dataIndex, : );
    thisEndPosSimple = MCSimplePosEnd( dataIndex, : );
    
    thisMeanEndPosHPOP = MC_HPOP_PosEndMean( dataIndex, : );
    thisEndPosHPOP = MC_HPOP_PosEnd( dataIndex, : );
    
    absMeanDeviationEndPosSimple( dataIndex ) = norm( rECIManouverEnd_target - thisMeanEndPosSimple' );
    absMeanDeviationEndPosHPOP( dataIndex ) = norm( rECIManouverEnd_target - thisMeanEndPosHPOP' );
    
    relEndPosSimple_chaser( dataIndex, : ) = relEndPosSimple_chaser( dataIndex, : ) - rECIManouverEnd_target';
    relEndPosHPOP_chaser( dataIndex, : ) = relEndPosHPOP_chaser( dataIndex, : ) - rECIManouverEnd_target';
    
    absDeviationEndPosSimple( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosSimple' );
    absDeviationEndPosHPOP( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosHPOP' );
    
end



%% Plots


% figure(1)
% hold on
% title('Norm of Error in Point of Rendezvous of Mean Value MC Simulations')
% plot( absMeanDeviationEndPosSimple )
% plot( absMeanDeviationEndPosHPOP )
% legend('Simple Model', 'HPOP')
% hold off
% 
% figure(2)
% hold on
% title('Error in ECI XY-Position of MC Simulations')
% plot( relEndPosSimple_chaser( : , 1 ), relEndPosSimple_chaser( :, 2 ), '*' )
% plot( relEndPosHPOP_chaser( : , 1 ), relEndPosHPOP_chaser( :, 2 ), '*' )
% legend('Simple Model', 'HPOP')
% xlabel('X')
% ylabel('Y')
% hold off
% 
% figure(3)
% hold on
% title('Error in ECI XZ-Position of MC Simulations')
% plot( relEndPosSimple_chaser( :, 1 ), relEndPosSimple_chaser( :, 3 ), '*' )
% plot( relEndPosHPOP_chaser( :, 1 ), relEndPosHPOP_chaser( :, 3 ), '*' )
% legend('Simple Model', 'HPOP')
% xlabel('X')
% ylabel('Z')
% hold off
% 
% figure(4)
% hold on
% title('Error in ECI XYZ-Position of MC Simulations')
% plot3( relEndPosSimple_chaser( :, 1 ), relEndPosSimple_chaser( :, 2 ), relEndPosSimple_chaser( :, 3 ), '*' )
% plot3( relEndPosHPOP_chaser( :, 1 ), relEndPosHPOP_chaser( :, 2 ), relEndPosHPOP_chaser( :, 3 ), '*' )
% legend('Simple Model', 'HPOP')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% hold off
% 
% figure(5)
% hold on
% title('Norm of Error in Point of Rendezvous of MC Simulations')
% plot( absDeviationEndPosSimple, '*' )
% plot( absDeviationEndPosHPOP, '*' )
% legend('Simple Model', 'HPOP')
% hold off
% 
% figure(6)
% hold on
% title('Time Delays of MC Simulations')
% plot( MCdeviationTimes, '*' )
% hold off
% 
% figure(7)
% hold on
% title('ECI Trajectories')
% plot3(rXECITrajectoryInitial_chaser, rYECITrajectoryInitial_chaser, rZECITrajectoryInitial_chaser)
% plot3( r0ECI_chaser(1), r0ECI_chaser(2), r0ECI_chaser(3), '*' )
% plot3( rECIExperimentStartSimple_chaser(1), rECIExperimentStartSimple_chaser(2), rECIExperimentStartSimple_chaser(3), '*' )
% plot3( rECIExperimentEndSimple_chaser(1), rECIExperimentEndSimple_chaser(2), rECIExperimentEndSimple_chaser(3), '*' )
% plot3( Eph(:, 2)./10^3, Eph(:, 3)./10^3, Eph(:, 4)./10^3)
% plot3(rXECITrajectoryNew_chaser, rYECITrajectoryNew_chaser, rZECITrajectoryNew_chaser)
% legend('Simple Model initial trajectory', 'Simple Model initial point', 'Simple Model maneuver start', 'Simple Model maneuver end', 'HPOP', 'Simple Model maneuver trajectory')
% hold off
% 
figure(8)
hold on
title('Norm of Error in Point of Rendezvous HPOP MC Simulation')
plot( absDeviationEndPosHPOP )
hold off

% figure(9)
% hold on
% title('Norm of Error in Point of Rendezvous of Mean Value HPOP MC Simulation')
% plot( absMeanDeviationEndPosHPOP )
% hold off

figure(10)
hold on
grid on
title('ECI HPOP Trajectories')
plot3( r0ECI_chaser(1), r0ECI_chaser(2), r0ECI_chaser(3), '*' )
plot3( rECIExperimentStartSimple_chaser(1), rECIExperimentStartSimple_chaser(2), rECIExperimentStartSimple_chaser(3), '*' )
plot3( rECIExperimentEndSimple_chaser(1), rECIExperimentEndSimple_chaser(2), rECIExperimentEndSimple_chaser(3), '*' )
for plotIndex = 1 : MCsampleNum
    plot3( MC_HPOP_ECI_X_Trajectories(plotIndex, :), MC_HPOP_ECI_Y_Trajectories(plotIndex, :), MC_HPOP_ECI_Z_Trajectories(plotIndex, :))
end
legend('Simple Model initial point', 'Simple Model maneuver start', 'Simple Model maneuver end')
hold off

figure(11)
hold on
grid on
title('ECI Simple Trajectories')
plot3( r0ECI_chaser(1), r0ECI_chaser(2), r0ECI_chaser(3), '*' )
plot3( rECIExperimentStartSimple_chaser(1), rECIExperimentStartSimple_chaser(2), rECIExperimentStartSimple_chaser(3), '*' )
plot3( rECIExperimentEndSimple_chaser(1), rECIExperimentEndSimple_chaser(2), rECIExperimentEndSimple_chaser(3), '*' )
for plotIndex = 1 : MCsampleNum
    plot3( MCSimpleECI_X_Trajectories(plotIndex, :), MCSimpleECI_Y_Trajectories(plotIndex, :), MCSimpleECI_Z_Trajectories(plotIndex, :))
end
legend('Simple Model initial point', 'Simple Model maneuver start', 'Simple Model maneuver end')
hold off

figure(12)
hold on
title('Norm of Error in Point of Rendezvous Simple Model MC Simulation')
plot( absDeviationEndPosSimple )
hold off
