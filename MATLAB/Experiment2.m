% Experiment 2 is based on Example 7.4 in Curtis
% Clear Workspace
clc
clear all
format long g
close all

%% Load Parameters

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC

run earthParameters; 
run satelliteParametersExperiment2

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0,...
                  'Thrust',0, 'stepCounter', 0, 'thrustStartTime', 0, 'thrustEndTime', 0, 'VelocityChange', 0, 'thrustAcceleration', 0);

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 2000;

run SAT_Const
run constants
load DE436Coeff.mat
PC = DE436Coeff; % JPL Planetary Ephemeris DECEMBER 14 1949 - DECEMBER 21 2149

run loadHPOPFileData

% epoch state 
year = 2002;
mon = 04;
day = 24;
hour = 00;
min = 00;
sec = 00;
Y0(1) = rTargetECI(1)*10^3;
Y0(2) = rTargetECI(2)*10^3;
Y0(3) = rTargetECI(3)*10^3;
Y0(4) = vTargetECI(1)*10^3;
Y0(5) = vTargetECI(2)*10^3;
Y0(6) = vTargetECI(3)*10^3;
AuxParam.area_solar = areaSolarChaser;
AuxParam.area_drag = areaDragChaser;
AuxParam.mass = massChaser;
AuxParam.Cr = CrChaser;
AuxParam.Cd = CdChaser;

% epoch
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
%Y0 = ECEF2ECI(Mjd_UTC, Y0); Already in ECI


%% Experiment Setup

maneuverEndTime = 0.25*3600 + 4380; %4875;
maneuverStartTime = 0.25*3600;
AuxParam.thrustDuration = 10;
Step   = 0.1;   % [s]
N_Step = maneuverEndTime*1/Step; % 

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

orbitType_chaser = "prograde";

numSamplePointsInitialTrajectorySimple = maneuverStartTime;
numSamplePointsFinalTrajectorySimple = maneuverEndTime - maneuverStartTime;

%% Monte Carlo Experiment Setup

MCsampleNum = 10;

meanDeviationTimeSetup = 0;
maxDeviationTimeSetup = 2;

MCdeviationTimes = meanDeviationTimeSetup - maxDeviationTimeSetup + ( 2 * maxDeviationTimeSetup * rand( MCsampleNum, 1 ) );

experimentStartTimeIdeal = maneuverStartTime;


%% Initial Orbit Determination

r0ECI_target = rTargetECI;
v0ECI_target = vTargetECI;

r0ECI_chaser = rChaserECI;
v0ECI_chaser = vChaserECI;

orbitPeriod_chaser = tChaser;

% A's position at experiment start
[ rECIExperimentStart_target, vECIExperimentStart_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, experimentStartTimeIdeal, anomalyErrorTolerance, anomalyMaxIterations );
[rXECITrajectoryTarget, rYECITrajectoryTarget, rZECITrajectoryTarget, vXECITrajectoryTarget, vYECITrajectoryTarget, vZECITrajectoryTarget, sampleTECITrajectoryTarget] = ECITrajectory( r0ECI_target, v0ECI_target, anomalyErrorTolerance, anomalyMaxIterations, maneuverEndTime, 1, N_Step+1, muEarth );

% A's position after ideal manouver
[ rECIManouverEnd_target, vECIManouverEnd_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverEndTime, anomalyErrorTolerance, anomalyMaxIterations );
 
% B's ideal position at experiment start
[ rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStartECI_chaser, deltaVEndECI_chaser, vIntersectOrbit_chaser ] = interceptOrbit( rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser, rECIManouverEnd_target, vECIManouverEnd_target, maneuverEndTime - maneuverStartTime, orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

[ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser );
deltaVStartLVLH_chaser = QmatECItoLVLH_chaser * deltaVStartECI_chaser;
deltaVEndLVLH_chaser = QmatECItoLVLH_chaser * deltaVEndECI_chaser;

%% Find optimal maneuver time

% 
% minManeuverTime = 60; % [s]
% maxManeuverTime = maneuverEndTime;
% deltaManeuverTime = 1;
% optimalManeuverTime = minManeuverTime;
% optimalManeuverDeltaV = 10^9;
% 
% deltaVManeuverTimes = zeros( 1 + ( maxManeuverTime - minManeuverTime ) / deltaManeuverTime, 1 );
% deltaVStartManeuverTimes = zeros( 1 + ( maxManeuverTime - minManeuverTime ) / deltaManeuverTime, 1 );
% deltaVEndManeuverTimes = zeros( 1 + ( maxManeuverTime - minManeuverTime ) / deltaManeuverTime, 1 );
% 
% 
% for thisManeuverTime = minManeuverTime : deltaManeuverTime : maxManeuverTime
%     
%     % A's position after ideal manouver
%     [ rTarget, vTarget ] = nextStateTimeStep( muEarth, rECIExperimentStart_target, vECIExperimentStart_target, thisManeuverTime, anomalyErrorTolerance, anomalyMaxIterations );
% 
%     % Required velocity change satellite B
%     [ deltaVStart, deltaVEnd ] = interceptOrbit( rECIExperimentStartIdeal_chaser, vECIExperimentStartIdeal_chaser, rTarget, vTarget, thisManeuverTime, orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
% 
%     thisManeuverDeltaV = norm( deltaVStart ) + norm( deltaVEnd );
%     
%     if norm( deltaVStart ) < optimalManeuverDeltaV
%         
%         optimalManeuverTime = thisManeuverTime;
%         optimalManeuverDeltaV = norm( deltaVStart );
%         
%     end
%         
%     deltaVManeuverTimes( 1 + (thisManeuverTime - minManeuverTime ) / deltaManeuverTime ) = thisManeuverDeltaV;
%     deltaVStartManeuverTimes( 1 + (thisManeuverTime - minManeuverTime ) / deltaManeuverTime ) = norm( deltaVStart );
%     deltaVEndManeuverTimes( 1 + (thisManeuverTime - minManeuverTime ) / deltaManeuverTime ) = norm( deltaVEnd );
%     
% end
% 
% figure(1)
% hold on
% grid on
% title('Start + End deltaV')
% plot( minManeuverTime : deltaManeuverTime : maxManeuverTime, deltaVManeuverTimes' )
% hold off
% 
% figure(2)
% hold on
% grid on
% title('Start deltaV')
% plot( minManeuverTime : deltaManeuverTime : maxManeuverTime, deltaVStartManeuverTimes' )
% hold off
% 
% figure(3)
% hold on
% grid on
% title('End deltaV')
% plot( minManeuverTime : deltaManeuverTime : maxManeuverTime, deltaVEndManeuverTimes' )
% hold off

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


%%%%%%%%%%%%%  HPOP MODEL
% propagation
AuxParam.thrustAcceleration = (deltaVExperimentStart_chaser*1000)./AuxParam.thrustDuration;
AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + MCdeviationTimes( 1 ))/86400);
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

    
    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    AuxParam.thrustAcceleration = (deltaVManeuverStart_chaser*1000)./AuxParam.thrustDuration;
    AuxParam.thrustStartTime = Mjd0 + ((maneuverStartTime + thisDeviationTime)/86400);
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

relXTrajectoriesHPOP = zeros( MCsampleNum, N_Step+1 );
relYTrajectoriesHPOP = zeros( MCsampleNum, N_Step+1 ); 
relZTrajectoriesHPOP = zeros( MCsampleNum, N_Step+1 );
relXTrajectoriesSimple = zeros( MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);
relYTrajectoriesSimple = zeros( MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);
relZTrajectoriesSimple = zeros( MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);

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

for MCIndex = 1 : MCsampleNum
    
    for HPOPIndex = 1 : N_Step + 1
        
        ECIPosTarget = [ rXECITrajectoryTarget( HPOPIndex ), rYECITrajectoryTarget( HPOPIndex ), rZECITrajectoryTarget( HPOPIndex )]';
        
        ECIVelTarget = [ vXECITrajectoryTarget( HPOPIndex ), vYECITrajectoryTarget( HPOPIndex ), vZECITrajectoryTarget( HPOPIndex )]';
        
        ECIPosChaser = [ MC_HPOP_ECI_X_Trajectories( MCIndex, HPOPIndex ), MC_HPOP_ECI_Y_Trajectories( MCIndex, HPOPIndex ), MC_HPOP_ECI_Z_Trajectories( MCIndex, HPOPIndex )]';
        
        ChaserRelativeToTarget = BPosRelativeToA( ECIPosTarget, ECIVelTarget, ECIPosChaser );
        
        relXTrajectoriesHPOP( MCIndex, HPOPIndex ) = ChaserRelativeToTarget(1);
        relYTrajectoriesHPOP( MCIndex, HPOPIndex ) = ChaserRelativeToTarget(2);
        relZTrajectoriesHPOP( MCIndex, HPOPIndex ) = ChaserRelativeToTarget(3);
        
    end
    
    for SimpleIndex = 1 : numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple
        
        ECIPosTarget = [ rXECITrajectoryTarget( SimpleIndex*10 ), rYECITrajectoryTarget( SimpleIndex*10 ), rZECITrajectoryTarget( SimpleIndex*10 )]';
        
        ECIVelTarget = [ vXECITrajectoryTarget( SimpleIndex*10 ), vYECITrajectoryTarget( SimpleIndex*10 ), vZECITrajectoryTarget( SimpleIndex*10 )]';
        
        ECIPosChaser = [ MCSimpleECI_X_Trajectories( MCIndex, SimpleIndex ), MCSimpleECI_Z_Trajectories( MCIndex, SimpleIndex ), MCSimpleECI_Z_Trajectories( MCIndex, SimpleIndex )]';
        
        ChaserRelativeToTarget = BPosRelativeToA( ECIPosTarget, ECIVelTarget, ECIPosChaser );
        
        relXTrajectoriesSimple( MCIndex, HPOPIndex ) = ChaserRelativeToTarget(1);
        relYTrajectoriesSimple( MCIndex, HPOPIndex ) = ChaserRelativeToTarget(2);
        relZTrajectoriesSimple( MCIndex, HPOPIndex ) = ChaserRelativeToTarget(3);
        
    end
    
end


%% Plots


figure(1)
hold on
title('Norm of Error in Point of Rendezvous of Mean Value MC Simulations')
plot( absMeanDeviationEndPosSimple )
plot( absMeanDeviationEndPosHPOP )
legend('Simple Model', 'HPOP')
hold off

figure(2)
hold on
title('Error in ECI XY-Position of MC Simulations')
plot( relEndPosSimple_chaser( : , 1 ), relEndPosSimple_chaser( :, 2 ), '*' )
plot( relEndPosHPOP_chaser( : , 1 ), relEndPosHPOP_chaser( :, 2 ), '*' )
legend('Simple Model', 'HPOP')
xlabel('X')
ylabel('Y')
hold off

figure(3)
hold on
title('Error in ECI XZ-Position of MC Simulations')
plot( relEndPosSimple_chaser( :, 1 ), relEndPosSimple_chaser( :, 3 ), '*' )
plot( relEndPosHPOP_chaser( :, 1 ), relEndPosHPOP_chaser( :, 3 ), '*' )
legend('Simple Model', 'HPOP')
xlabel('X')
ylabel('Z')
hold off

figure(4)
hold on
title('Error in ECI XYZ-Position of MC Simulations')
plot3( relEndPosSimple_chaser( :, 1 ), relEndPosSimple_chaser( :, 2 ), relEndPosSimple_chaser( :, 3 ), '*' )
plot3( relEndPosHPOP_chaser( :, 1 ), relEndPosHPOP_chaser( :, 2 ), relEndPosHPOP_chaser( :, 3 ), '*' )
legend('Simple Model', 'HPOP')
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off

figure(5)
hold on
title('Norm of Error in Point of Rendezvous of MC Simulations')
plot( absDeviationEndPosSimple, '*' )
plot( absDeviationEndPosHPOP, '*' )
legend('Simple Model', 'HPOP')
hold off

figure(6)
hold on
title('Time Delays of MC Simulations')
plot( MCdeviationTimes, '*' )
hold off

figure(7)
hold on
title('ECI Trajectories')
plot3(rXECITrajectoryInitial_chaser, rYECITrajectoryInitial_chaser, rZECITrajectoryInitial_chaser)
plot3( r0ECI_chaser(1), r0ECI_chaser(2), r0ECI_chaser(3), '*' )
plot3( rECIExperimentStartSimple_chaser(1), rECIExperimentStartSimple_chaser(2), rECIExperimentStartSimple_chaser(3), '*' )
plot3( rECIExperimentEndSimple_chaser(1), rECIExperimentEndSimple_chaser(2), rECIExperimentEndSimple_chaser(3), '*' )
plot3( Eph(:, 2)./10^3, Eph(:, 3)./10^3, Eph(:, 4)./10^3)
plot3(rXECITrajectoryNew_chaser, rYECITrajectoryNew_chaser, rZECITrajectoryNew_chaser)
legend('Simple Model initial trajectory', 'Simple Model initial point', 'Simple Model maneuver start', 'Simple Model maneuver end', 'HPOP', 'Simple Model maneuver trajectory')
hold off

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

figure(13)
hold on
grid on
title('ECI Simple and HPOP Trajectories')
plot3( r0ECI_chaser(1), r0ECI_chaser(2), r0ECI_chaser(3), '*' )
plot3( rECIExperimentStartSimple_chaser(1), rECIExperimentStartSimple_chaser(2), rECIExperimentStartSimple_chaser(3), '*' )
plot3( rECIExperimentEndSimple_chaser(1), rECIExperimentEndSimple_chaser(2), rECIExperimentEndSimple_chaser(3), '*' )
for plotIndex = 1 : MCsampleNum
    plot3( MCSimpleECI_X_Trajectories(plotIndex, :), MCSimpleECI_Y_Trajectories(plotIndex, :), MCSimpleECI_Z_Trajectories(plotIndex, :), 'r')
end
for plotIndex = 1 : MCsampleNum
    plot3( MC_HPOP_ECI_X_Trajectories(plotIndex, :), MC_HPOP_ECI_Y_Trajectories(plotIndex, :), MC_HPOP_ECI_Z_Trajectories(plotIndex, :), 'b')
end
legend('Simple Model initial point', 'Simple Model maneuver start', 'Simple Model maneuver end')
hold off


figure(14)
hold on
grid on
title('Relative Trajectories MC Simple')
for plotIndex = 1 : MCsampleNum
    plot3( relXTrajectoriesSimple(plotIndex, :), relYTrajectoriesSimple(plotIndex, :), relZTrajectoriesSimple(plotIndex, :))
end
plot3( relXTrajectoriesSimple(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple), relYTrajectoriesSimple(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple), relZTrajectoriesSimple(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple), '*')
hold off

figure(15)
hold on
grid on
title('Relative Trajectories MC HPOP')
for plotIndex = 1 : MCsampleNum
    plot3( relXTrajectoriesHPOP(plotIndex, :), relYTrajectoriesHPOP(plotIndex, :), relZTrajectoriesHPOP(plotIndex, :))
end
plot3( relXTrajectoriesHPOP(MCsampleNum, N_Step + 1), relYTrajectoriesHPOP(MCsampleNum, N_Step + 1), relZTrajectoriesHPOP(MCsampleNum, N_Step + 1), '*')
hold off

figure(16)
hold on
grid on
title('Relative Trajectories MC HPOP and Simple')
for plotIndex = 1 : MCsampleNum
    plot3( relXTrajectoriesHPOP(plotIndex, :), relYTrajectoriesHPOP(plotIndex, :), relZTrajectoriesHPOP(plotIndex, :), 'r')
end
for plotIndex = 1 : MCsampleNum
    plot3( relXTrajectoriesSimple(plotIndex, :), relYTrajectoriesSimple(plotIndex, :), relZTrajectoriesSimple(plotIndex, :), 'b')
end
plot3( relXTrajectoriesSimple(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple), relYTrajectoriesSimple(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple), relZTrajectoriesSimple(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple), '*')
plot3( relXTrajectoriesHPOP(MCsampleNum, N_Step + 1), relYTrajectoriesHPOP(MCsampleNum, N_Step + 1), relZTrajectoriesHPOP(MCsampleNum, N_Step + 1), '*')
hold off
