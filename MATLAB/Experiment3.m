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
Step   = 0.1;   % [s]
N_Step = maneuverEndTime*1/Step; %             
              
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
%Cnm = Cnm(1:AuxParam.n+1,1:AuxParam.n+1);
%Snm = Snm(1:AuxParam.n+1,1:AuxParam.n+1);

orbitType_chaser = "retrograde";

numSamplePointsInitialTrajectorySimple = 200;
numSamplePointsFinalTrajectorySimple = 1000;

%% Monte Carlo Experiment Setup

MCsampleNum = 10;

meanDeviationTimeSetup = 0;
maxDeviationTimeSetup = 2;

MCdeviationTimes = meanDeviationTimeSetup - maxDeviationTimeSetup + ( 2 * maxDeviationTimeSetup * rand( MCsampleNum, 1 ) );

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


%% Monte Carlo Experimetns


%%%%%%%%%%%%%  SIMPLE MODEL


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


% (Debug) B's initial trajectory
[rXECITrajectoryInitial_chaser, rYECITrajectoryInitial_chaser, rZECITrajectoryInitial_chaser, vXECITrajectoryInitial_chaser, vYECITrajectoryInitial_chaser, vZECITrajectoryInitial_chaser, sampleTECITrajectoryInitial_chaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartTime, 1, numSamplePointsInitialTrajectorySimple, muEarth );

% B's position at first iteration start
[ rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + MCdeviationTimes( 1 ), anomalyErrorTolerance, anomalyMaxIterations );

[rXECITrajectoryExperimentStart_chaser, rYECITrajectoryExperimentStart_chaser, rZECITrajectoryExperimentStart_chaser, vXECITrajectoryExperimentStart_chaser, vYECITrajectoryExperimentStart_chaser, vZECITrajectoryExperimentStart_chaser, sampleTECITrajectoryExperimentStart_chaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartTime + MCdeviationTimes( 1 ), 1, numSamplePointsInitialTrajectorySimple, muEarth );
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
[ rECIExperimentEndSimple_chaser, vECIExperimentEndSimple_chaser ] = nextStateTimeStep( muEarth, rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser + deltaVExperimentStart_chaser, maneuverEndTime - maneuverStartTime - MCdeviationTimes( 1 ), anomalyErrorTolerance, anomalyMaxIterations );

% (Debug) B's new trajectory
[rXECITrajectoryNew_chaser, rYECITrajectoryNew_chaser, rZECITrajectoryNew_chaser, vXECITrajectoryNew_chaser, vYECITrajectoryNew_chaser, vZECITrajectoryNew_chaser, sampleTECITrajectoryNew_chaser] = ECITrajectory( rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser + deltaVExperimentStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverEndTime - maneuverStartTime, 1, numSamplePointsFinalTrajectorySimple, muEarth );

[rXECITrajectoryExperimentEnd_chaser, rYECITrajectoryExperimentEnd_chaser, rZECITrajectoryExperimentEnd_chaser, vXECITrajectoryExperimentEnd_chaser, vYECITrajectoryExperimentEnd_chaser, vZECITrajectoryExperimentEnd_chaser, sampleTECITrajectoryExperimentEnd_chaser] = ECITrajectory( rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser + deltaVExperimentStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverEndTime - maneuverStartTime - MCdeviationTimes( 1 ), 1, numSamplePointsFinalTrajectorySimple, muEarth );
MCSimpleECI_X_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rXECITrajectoryExperimentEnd_chaser;
MCSimpleECI_Y_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rYECITrajectoryExperimentEnd_chaser;
MCSimpleECI_Z_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rZECITrajectoryExperimentEnd_chaser;

MCSimplePosEnd( 1, : ) = rECIExperimentEndSimple_chaser';
MCSimpleVelEnd( 1, : ) = vECIExperimentEndSimple_chaser';

MCSimplePosEndMean( 1, : ) = MCSimplePosEnd( 1, : );
MCSimpleVelEndMean( 1, : ) = MCSimpleVelEnd( 1, : );



for experimentIndex = 2 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );
  
    %%%%%%%%%%%%%  SIMPLE MODEL
    % B's position at maneuver start
    [ rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    MCSimplePosStart( experimentIndex, : ) = rECIManeuverStartSimple_chaser';
    MCSimpleVelStart( experimentIndex, : ) = vECIManeuverStartSimple_chaser';
    
    [rXECITrajectoryManeuverStart_chaser, rYECITrajectoryManeuverStart_chaser, rZECITrajectoryManeuverStart_chaser, vXECITrajectoryManeuverStart_chaser, vYECITrajectoryManeuverStart_chaser, vZECITrajectoryManeuverStart_chaser, sampleTECITrajectoryManeuverStart_chaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartTime + MCdeviationTimes( experimentIndex ), 1, numSamplePointsInitialTrajectorySimple, muEarth );
    MCSimpleECI_X_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rXECITrajectoryManeuverStart_chaser;
    MCSimpleECI_Y_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rYECITrajectoryManeuverStart_chaser;
    MCSimpleECI_Z_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rZECITrajectoryManeuverStart_chaser;

    MCSimplePosStartMean( experimentIndex, : ) = mean( MCSimplePosStart( 1 : experimentIndex, : ) );
    
    [ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser );
    QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
    deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;

    % B's position at experiment end
    [ rECIExperimentEndSimple_chaser, vECIExperimentEndSimple_chaser ] = nextStateTimeStep( muEarth, rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser + deltaVManeuverStart_chaser, maneuverEndTime - maneuverStartTime - thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    MCSimplePosEnd( experimentIndex, : ) = rECIExperimentEndSimple_chaser';
    MCSimpleVelEnd( experimentIndex, : ) = vECIExperimentEndSimple_chaser';
    
    MCSimplePosEndMean( experimentIndex, : ) = mean( MCSimplePosEnd( 1 : experimentIndex, : ) );
    MCSimpleVelEndMean( experimentIndex, : ) = mean( MCSimpleVelEnd( 1 : experimentIndex, : ) );
    
    [rXECITrajectoryManeuverEnd_chaser, rYECITrajectoryManeuverEnd_chaser, rZECITrajectoryManeuverEnd_chaser, vXECITrajectoryManeuverEnd_chaser, vYECITrajectoryManeuverEnd_chaser, vZECITrajectoryManeuverEnd_chaser, sampleTECITrajectoryManeuverEnd_chaser] = ECITrajectory( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser + deltaVManeuverStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverEndTime - maneuverStartTime - MCdeviationTimes( experimentIndex ), 1, numSamplePointsFinalTrajectorySimple, muEarth );
    MCSimpleECI_X_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rXECITrajectoryManeuverEnd_chaser;
    MCSimpleECI_Y_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rYECITrajectoryManeuverEnd_chaser;
    MCSimpleECI_Z_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rZECITrajectoryManeuverEnd_chaser;
    
end


%%%%%%%%%%%%%  HPOP MODEL


%%%%%%%  Experiment 1 Only Thrust effect d = 10s


AuxParam.thrustDuration = 10;
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


MC_1_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_1_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_1_HPOP_PosEndMean = zeros(MCsampleNum, 3);
MC_1_HPOP_PosStartMean = zeros(MCsampleNum, 3);

MC_1_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_1_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_1_HPOP_VelEndMean = zeros(MCsampleNum, 3);
MC_1_HPOP_VelStartMean = zeros(MCsampleNum, 3);

MC_1_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_1_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_1_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);


% propagation
AuxParam.thrustAcceleration = (deltaVExperimentStart_chaser*1000)./AuxParam.thrustDuration;
AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + MCdeviationTimes( 1 ))/86400);
[Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
MC_1_HPOP_PosEnd( 1, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
MC_1_HPOP_VelEnd( 1, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

MC_1_HPOP_PosEndMean( 1, : ) = MC_1_HPOP_PosEnd( 1, : );
MC_1_HPOP_VelEndMean( 1, : ) = MC_1_HPOP_VelEnd( 1, : );

MC_1_HPOP_ECI_X_Trajectories(1, :) = Eph(:, 2)'./10^3;
MC_1_HPOP_ECI_Y_Trajectories(1, :) = Eph(:, 3)'./10^3;
MC_1_HPOP_ECI_Z_Trajectories(1, :) = Eph(:, 4)'./10^3;


for experimentIndex = 2 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );
    
    % B's position at maneuver start
    [ rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    [ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser );
    QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
    deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;
    
    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    AuxParam.thrustAcceleration = (deltaVManeuverStart_chaser*1000)./AuxParam.thrustDuration;
    AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + thisDeviationTime)/86400);
    [Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
    MC_1_HPOP_PosEnd( experimentIndex, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
    MC_1_HPOP_VelEnd( experimentIndex, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

    MC_1_HPOP_PosEndMean( experimentIndex, : ) = mean(MC_1_HPOP_PosEnd( 1 : experimentIndex, : ));
    MC_1_HPOP_VelEndMean( experimentIndex, : ) = mean(MC_1_HPOP_VelEnd( 1 : experimentIndex, : ));
    
    MC_1_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph(:, 2)'./10^3;
    MC_1_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph(:, 3)'./10^3;
    MC_1_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph(:, 4)'./10^3;
    
end



%%%%%%%  Experiment 1 Only Thrust effect d = 20s


AuxParam.thrustDuration = 20;
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


MC_2_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_2_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_2_HPOP_PosEndMean = zeros(MCsampleNum, 3);
MC_2_HPOP_PosStartMean = zeros(MCsampleNum, 3);

MC_2_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_2_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_2_HPOP_VelEndMean = zeros(MCsampleNum, 3);
MC_2_HPOP_VelStartMean = zeros(MCsampleNum, 3);

MC_2_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_2_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_2_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);


% propagation
AuxParam.thrustAcceleration = (deltaVExperimentStart_chaser*1000)./AuxParam.thrustDuration;
AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + MCdeviationTimes( 1 ))/86400);
[Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
MC_2_HPOP_PosEnd( 1, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
MC_2_HPOP_VelEnd( 1, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

MC_2_HPOP_PosEndMean( 1, : ) = MC_2_HPOP_PosEnd( 1, : );
MC_2_HPOP_VelEndMean( 1, : ) = MC_2_HPOP_VelEnd( 1, : );

MC_2_HPOP_ECI_X_Trajectories(1, :) = Eph(:, 2)'./10^3;
MC_2_HPOP_ECI_Y_Trajectories(1, :) = Eph(:, 3)'./10^3;
MC_2_HPOP_ECI_Z_Trajectories(1, :) = Eph(:, 4)'./10^3;


for experimentIndex = 2 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );
    
    % B's position at maneuver start
    [ rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    [ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser );
    QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
    deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;
    
    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    AuxParam.thrustAcceleration = (deltaVManeuverStart_chaser*1000)./AuxParam.thrustDuration;
    AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + thisDeviationTime)/86400);
    [Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
    MC_2_HPOP_PosEnd( experimentIndex, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
    MC_2_HPOP_VelEnd( experimentIndex, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

    MC_2_HPOP_PosEndMean( experimentIndex, : ) = mean(MC_2_HPOP_PosEnd( 1 : experimentIndex, : ));
    MC_2_HPOP_VelEndMean( experimentIndex, : ) = mean(MC_2_HPOP_VelEnd( 1 : experimentIndex, : ));
    
    MC_2_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph(:, 2)'./10^3;
    MC_2_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph(:, 3)'./10^3;
    MC_2_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph(:, 4)'./10^3;
    
end



%%%%%%%  Experiment 3 Thrust effect d = 10s, only J2 effects


AuxParam.thrustDuration = 10;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 2;
AuxParam.m       = 2;
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


MC_3_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_3_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_3_HPOP_PosEndMean = zeros(MCsampleNum, 3);
MC_3_HPOP_PosStartMean = zeros(MCsampleNum, 3);

MC_3_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_3_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_3_HPOP_VelEndMean = zeros(MCsampleNum, 3);
MC_3_HPOP_VelStartMean = zeros(MCsampleNum, 3);

MC_3_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_3_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_3_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);


% propagation
AuxParam.thrustAcceleration = (deltaVExperimentStart_chaser*1000)./AuxParam.thrustDuration;
AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + MCdeviationTimes( 1 ))/86400);
[Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
MC_3_HPOP_PosEnd( 1, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
MC_3_HPOP_VelEnd( 1, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

MC_3_HPOP_PosEndMean( 1, : ) = MC_3_HPOP_PosEnd( 1, : );
MC_3_HPOP_VelEndMean( 1, : ) = MC_3_HPOP_VelEnd( 1, : );

MC_3_HPOP_ECI_X_Trajectories(1, :) = Eph(:, 2)'./10^3;
MC_3_HPOP_ECI_Y_Trajectories(1, :) = Eph(:, 3)'./10^3;
MC_3_HPOP_ECI_Z_Trajectories(1, :) = Eph(:, 4)'./10^3;


for experimentIndex = 2 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );
    
    % B's position at maneuver start
    [ rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    [ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser );
    QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
    deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;
    
    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    AuxParam.thrustAcceleration = (deltaVManeuverStart_chaser*1000)./AuxParam.thrustDuration;
    AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + thisDeviationTime)/86400);
    [Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
    MC_3_HPOP_PosEnd( experimentIndex, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
    MC_3_HPOP_VelEnd( experimentIndex, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

    MC_3_HPOP_PosEndMean( experimentIndex, : ) = mean(MC_3_HPOP_PosEnd( 1 : experimentIndex, : ));
    MC_3_HPOP_VelEndMean( experimentIndex, : ) = mean(MC_3_HPOP_VelEnd( 1 : experimentIndex, : ));
    
    MC_3_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph(:, 2)'./10^3;
    MC_3_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph(:, 3)'./10^3;
    MC_3_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph(:, 4)'./10^3;
    
end




%%%%%%%  Experiment 4 Thrust effect d = 10s, all J2...40 effects, drag,
%%%%%%%  radiation pressure, sun, moon, planets, ocean tides


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
AuxParam.VelocityChange = 0;  


MC_4_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_4_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_4_HPOP_PosEndMean = zeros(MCsampleNum, 3);
MC_4_HPOP_PosStartMean = zeros(MCsampleNum, 3);

MC_4_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_4_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_4_HPOP_VelEndMean = zeros(MCsampleNum, 3);
MC_4_HPOP_VelStartMean = zeros(MCsampleNum, 3);

MC_4_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_4_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_4_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);


% propagation
AuxParam.thrustAcceleration = (deltaVExperimentStart_chaser*1000)./AuxParam.thrustDuration;
AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + MCdeviationTimes( 1 ))/86400);
[Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
MC_4_HPOP_PosEnd( 1, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
MC_4_HPOP_VelEnd( 1, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

MC_4_HPOP_PosEndMean( 1, : ) = MC_4_HPOP_PosEnd( 1, : );
MC_4_HPOP_VelEndMean( 1, : ) = MC_4_HPOP_VelEnd( 1, : );

MC_4_HPOP_ECI_X_Trajectories(1, :) = Eph(:, 2)'./10^3;
MC_4_HPOP_ECI_Y_Trajectories(1, :) = Eph(:, 3)'./10^3;
MC_4_HPOP_ECI_Z_Trajectories(1, :) = Eph(:, 4)'./10^3;


for experimentIndex = 2 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );
    
    % B's position at maneuver start
    [ rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    [ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser );
    QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
    deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;
    
    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    AuxParam.thrustAcceleration = (deltaVManeuverStart_chaser*1000)./AuxParam.thrustDuration;
    AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + thisDeviationTime)/86400);
    [Eph, stats] = ephemeris_Experiment1(Y0, N_Step, Step);
    MC_4_HPOP_PosEnd( experimentIndex, : ) = [Eph(N_Step+1, 2), Eph(N_Step+1, 3), Eph(N_Step+1, 4)]./10^3;
    MC_4_HPOP_VelEnd( experimentIndex, : ) = [Eph(N_Step+1, 5), Eph(N_Step+1, 6), Eph(N_Step+1, 7)]./10^3;

    MC_4_HPOP_PosEndMean( experimentIndex, : ) = mean(MC_4_HPOP_PosEnd( 1 : experimentIndex, : ));
    MC_4_HPOP_VelEndMean( experimentIndex, : ) = mean(MC_4_HPOP_VelEnd( 1 : experimentIndex, : ));
    
    MC_4_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph(:, 2)'./10^3;
    MC_4_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph(:, 3)'./10^3;
    MC_4_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph(:, 4)'./10^3;
    
end



%% Statistical Analysis

absMeanDeviationEndPosSimple = zeros( MCsampleNum, 1 );
absDeviationEndPosSimple = zeros( MCsampleNum, 1 );
relEndPosSimple_chaser = MCSimplePosEnd;

absMeanDeviationEndPosHPOP_1 = zeros( MCsampleNum, 1 );
absDeviationEndPosHPOP_1 = zeros( MCsampleNum, 1 );
relEndPosHPOP_1_chaser = MC_1_HPOP_PosEnd;

absMeanDeviationEndPosHPOP_2 = zeros( MCsampleNum, 1 );
absDeviationEndPosHPOP_2 = zeros( MCsampleNum, 1 );
relEndPosHPOP_2_chaser = MC_2_HPOP_PosEnd;

absMeanDeviationEndPosHPOP_3 = zeros( MCsampleNum, 1 );
absDeviationEndPosHPOP_3 = zeros( MCsampleNum, 1 );
relEndPosHPOP_3_chaser = MC_3_HPOP_PosEnd;

absMeanDeviationEndPosHPOP_4 = zeros( MCsampleNum, 1 );
absDeviationEndPosHPOP_4 = zeros( MCsampleNum, 1 );
relEndPosHPOP_4_chaser = MC_4_HPOP_PosEnd;

for dataIndex = 1 : MCsampleNum
    
    thisMeanEndPosSimple = MCSimplePosEndMean( dataIndex, : );
    thisEndPosSimple = MCSimplePosEnd( dataIndex, : );
    
    thisMeanEndPosHPOP_1 = MC_1_HPOP_PosEndMean( dataIndex, : );
    thisEndPosHPOP_1 = MC_1_HPOP_PosEnd( dataIndex, : );
    
    thisMeanEndPosHPOP_2 = MC_2_HPOP_PosEndMean( dataIndex, : );
    thisEndPosHPOP_2 = MC_2_HPOP_PosEnd( dataIndex, : );
    
    thisMeanEndPosHPOP_3 = MC_3_HPOP_PosEndMean( dataIndex, : );
    thisEndPosHPOP_3 = MC_3_HPOP_PosEnd( dataIndex, : );
    
    thisMeanEndPosHPOP_4 = MC_4_HPOP_PosEndMean( dataIndex, : );
    thisEndPosHPOP_4 = MC_4_HPOP_PosEnd( dataIndex, : );
    
    absMeanDeviationEndPosSimple( dataIndex ) = norm( rECIManouverEnd_target - thisMeanEndPosSimple' );
    absMeanDeviationEndPosHPOP_1( dataIndex ) = norm( rECIManouverEnd_target - thisMeanEndPosHPOP_1' );
    absMeanDeviationEndPosHPOP_2( dataIndex ) = norm( rECIManouverEnd_target - thisMeanEndPosHPOP_2' );
    absMeanDeviationEndPosHPOP_3( dataIndex ) = norm( rECIManouverEnd_target - thisMeanEndPosHPOP_3' );
    absMeanDeviationEndPosHPOP_4( dataIndex ) = norm( rECIManouverEnd_target - thisMeanEndPosHPOP_4' );
    
    relEndPosSimple_chaser( dataIndex, : ) = relEndPosSimple_chaser( dataIndex, : ) - rECIManouverEnd_target';
    relEndPosHPOP_1_chaser( dataIndex, : ) = relEndPosHPOP_1_chaser( dataIndex, : ) - rECIManouverEnd_target';
    relEndPosHPOP_2_chaser( dataIndex, : ) = relEndPosHPOP_2_chaser( dataIndex, : ) - rECIManouverEnd_target';
    relEndPosHPOP_3_chaser( dataIndex, : ) = relEndPosHPOP_3_chaser( dataIndex, : ) - rECIManouverEnd_target';
    relEndPosHPOP_4_chaser( dataIndex, : ) = relEndPosHPOP_4_chaser( dataIndex, : ) - rECIManouverEnd_target';
   
    absDeviationEndPosSimple( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosSimple' );
    absDeviationEndPosHPOP_1( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosHPOP_1' );
    absDeviationEndPosHPOP_2( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosHPOP_2' );
    absDeviationEndPosHPOP_3( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosHPOP_3' );
    absDeviationEndPosHPOP_4( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosHPOP_4' );
    
end



%% Plots



figure(1)
hold on
grid on
title('Error in ECI XYZ-Position of MC Simulations')
plot3( relEndPosSimple_chaser( :, 1 ), relEndPosSimple_chaser( :, 2 ), relEndPosSimple_chaser( :, 3 ), '*' )
plot3( relEndPosHPOP_1_chaser( :, 1 ), relEndPosHPOP_1_chaser( :, 2 ), relEndPosHPOP_1_chaser( :, 3 ), '*' )
plot3( relEndPosHPOP_2_chaser( :, 1 ), relEndPosHPOP_2_chaser( :, 2 ), relEndPosHPOP_2_chaser( :, 3 ), '*' )
plot3( relEndPosHPOP_3_chaser( :, 1 ), relEndPosHPOP_3_chaser( :, 2 ), relEndPosHPOP_3_chaser( :, 3 ), '*' )
plot3( relEndPosHPOP_4_chaser( :, 1 ), relEndPosHPOP_4_chaser( :, 2 ), relEndPosHPOP_4_chaser( :, 3 ), '*' )
legend('Simple Model', 'HPOP 1', 'HPOP 2', 'HPOP 3', 'HPOP 4')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
hold off


figure(2)
hold on
grid on
title('Norm of Error in Point of Rendezvous of MC Simulations')
plot( absDeviationEndPosSimple, '*' )
plot( absDeviationEndPosHPOP_1, '*' )
plot( absDeviationEndPosHPOP_2, '*' )
plot( absDeviationEndPosHPOP_3, '*' )
plot( absDeviationEndPosHPOP_4, '*' )
legend('Simple Model', 'HPOP 1', 'HPOP 2', 'HPOP 3', 'HPOP 4')
xlabel('Sample')
ylabel('Distance [km]')
hold off


figure(3)
hold on
grid on
title('ECI Simple Trajectories')
plot3( r0ECI_chaser(1), r0ECI_chaser(2), r0ECI_chaser(3), '*m' )
plot3( rECIExperimentStartSimple_chaser(1), rECIExperimentStartSimple_chaser(2), rECIExperimentStartSimple_chaser(3), '*k' )
plot3( rECIExperimentEndSimple_chaser(1), rECIExperimentEndSimple_chaser(2), rECIExperimentEndSimple_chaser(3), '*k' )
for plotIndex = 1 : MCsampleNum
    plot3( MCSimpleECI_X_Trajectories(plotIndex, :), MCSimpleECI_Y_Trajectories(plotIndex, :), MCSimpleECI_Z_Trajectories(plotIndex, :), 'b')
    plot3( MC_1_HPOP_ECI_X_Trajectories(plotIndex, :), MC_1_HPOP_ECI_Y_Trajectories(plotIndex, :), MC_1_HPOP_ECI_Z_Trajectories(plotIndex, :), 'r')
    plot3( MC_2_HPOP_ECI_X_Trajectories(plotIndex, :), MC_2_HPOP_ECI_Y_Trajectories(plotIndex, :), MC_2_HPOP_ECI_Z_Trajectories(plotIndex, :), 'g')
    plot3( MC_3_HPOP_ECI_X_Trajectories(plotIndex, :), MC_3_HPOP_ECI_Y_Trajectories(plotIndex, :), MC_3_HPOP_ECI_Z_Trajectories(plotIndex, :), 'y')
    plot3( MC_4_HPOP_ECI_X_Trajectories(plotIndex, :), MC_4_HPOP_ECI_Y_Trajectories(plotIndex, :), MC_4_HPOP_ECI_Z_Trajectories(plotIndex, :), 'c')
end
legend('Simple Model initial point', 'Simple Model maneuver start', 'Simple Model maneuver end')
hold off


