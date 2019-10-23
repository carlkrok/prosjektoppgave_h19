%% Clear Workspace

clear all
close all

%% Load Parameters

run earthParameters; 
run satelliteParameters;

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 2000;

orbitType = "prograde";
orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );

maneuverTime = 2500;


%% Monte Carlo Experiment Setup

MCsampleNum = 1000;

meanDeviationTimeSetup = 0;
maxDeviationTimeSetup = 10;

MCdeviationTimes = meanDeviationTimeSetup - maxDeviationTimeSetup + ( 2 * maxDeviationTimeSetup * rand( MCsampleNum, 1 ) );

experimentStartTimeIdeal = 10 + 2 * maxDeviationTimeSetup;


%% Initial Orbit Determination

r0PQW_A = positionVectorPQW( muEarth, hNorm_A, e_A, T_A );
v0PQW_A = velocityVectorPQW( muEarth, hNorm_A, e_A, T_A );
QmatPQWtoECI_A = transformPQWtoECI( i_A, O_A, w_A );
r0ECI_A = QmatPQWtoECI_A * r0PQW_A;
v0ECI_A = QmatPQWtoECI_A * v0PQW_A;

r0PQW_B = positionVectorPQW( muEarth, hNorm_B, e_B, T_B );
v0PQW_B = velocityVectorPQW( muEarth, hNorm_B, e_B, T_B );
QmatPQWtoECI_B = transformPQWtoECI( i_B, O_B, w_B );
r0ECI_B = QmatPQWtoECI_B * r0PQW_B;
v0ECI_B = QmatPQWtoECI_B * v0PQW_B;

% A's position at experiment start
[ rECIExperimentStart_A, vECIExperimentStart_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, experimentStartTimeIdeal, anomalyErrorTolerance, anomalyMaxIterations );
    
% A's position after ideal manouver
[ rECIManouverEnd_A, vECIManouverEnd_A ] = nextStateTimeStep( muEarth, rECIExperimentStart_A, vECIExperimentStart_A, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );
 
% B's ideal position at experiment start
[ rECIExperimentStartIdeal_B, vECIExperimentStartIdeal_B ] = nextStateTimeStep( muEarth, r0ECI_B, v0ECI_B, experimentStartTimeIdeal, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStartECI_B, deltaVEndECI_B, vIntersectOrbit_B ] = interceptOrbit( rECIExperimentStartIdeal_B, vECIExperimentStartIdeal_B, rECIManouverEnd_A, vECIManouverEnd_A, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

[ QmatECItoLVLH_B ] = ECIToLVLH( rECIExperimentStartIdeal_B, vECIExperimentStartIdeal_B );
deltaVStartLVLH_B = QmatECItoLVLH_B * deltaVStartECI_B;
deltaVEndLVLH_B = QmatECItoLVLH_B * deltaVEndECI_B;


%% Monte Carlo Experimetns

MCposEnd = zeros( MCsampleNum, 3);
MCposStart = zeros( MCsampleNum, 3);
MCposEndMean = zeros(MCsampleNum, 3);
MCposStartMean = zeros(MCsampleNum, 3);

MCvelEnd = zeros(MCsampleNum, 3);
MCvelStart = zeros(MCsampleNum, 3);
MCvelEndMean = zeros(MCsampleNum, 3);
MCvelStartMean = zeros(MCsampleNum, 3);

% B's position at experiment start
[ rECIExperimentStart_B, vECIExperimentStart_B ] = nextStateTimeStep( muEarth, r0ECI_B, v0ECI_B, experimentStartTimeIdeal + MCdeviationTimes( 1 ), anomalyErrorTolerance, anomalyMaxIterations );

MCposStart( 1, : ) = rECIExperimentStart_B';
MCvelStart( 1, : ) = vECIExperimentStart_B';

MCposStartMean( 1, : ) = MCposStart( 1, : );

[ QmatECItoLVLH_B ] = ECIToLVLH( rECIExperimentStart_B, vECIExperimentStart_B );
QmatLVLHtoECI_B = QmatECItoLVLH_B';
deltaVExperimentStart_B = QmatLVLHtoECI_B * deltaVStartLVLH_B;

% B's position at experiment end
[ rECIExperimentEnd_B, vECIExperimentEnd_B ] = nextStateTimeStep( muEarth, rECIExperimentStart_B, vECIExperimentStart_B + deltaVExperimentStart_B, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );

MCposEnd( 1, : ) = rECIExperimentEnd_B';
MCvelEnd( 1, : ) = vECIExperimentEnd_B';

MCposEndMean( 1, : ) = MCposEnd( 1, : );
MCvelEndMean( 1, : ) = MCvelEnd( 1, : );

for experimentIndex = 2 : MCsampleNum
    
    thisDeviation = MCdeviationTimes( experimentIndex );
    
    % B's position at maneuver start
    [ rECIManeuverStart_B, vECIManeuverStart_B ] = nextStateTimeStep( muEarth, r0ECI_B, v0ECI_B, experimentStartTimeIdeal + thisDeviation, anomalyErrorTolerance, anomalyMaxIterations );

    MCposStart( experimentIndex, : ) = rECIManeuverStart_B';
    MCvelStart( experimentIndex, : ) = vECIManeuverStart_B';
    
    MCposStartMean( experimentIndex, : ) = mean( MCposStart( 1 : experimentIndex, : ) );
    
    [ QmatECItoLVLH_B ] = ECIToLVLH( rECIManeuverStart_B, vECIManeuverStart_B );
    QmatLVLHtoECI_B = QmatECItoLVLH_B';
    deltaVManeuverStart_B = QmatLVLHtoECI_B * deltaVStartLVLH_B;

    % B's position at experiment end
    [ rECIExperimentEnd_B, vECIExperimentEnd_B ] = nextStateTimeStep( muEarth, rECIManeuverStart_B, vECIManeuverStart_B + deltaVManeuverStart_B, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );

    MCposEnd( experimentIndex, : ) = rECIExperimentEnd_B';
    MCvelEnd( experimentIndex, : ) = vECIExperimentEnd_B';
    
    MCposEndMean( experimentIndex, : ) = mean( MCposEnd( 1 : experimentIndex, : ) );
    MCvelEndMean( experimentIndex, : ) = mean( MCvelEnd( 1 : experimentIndex, : ) );
    
end


%% Statistical Analysis

absDeviationEndPos = zeros( MCsampleNum, 1 );


for dataIndex = 1 : MCsampleNum
    
    thisMeanEndPos = MCposEndMean( dataIndex, : );
    
    absDeviationEndPos( dataIndex ) = norm( rECIManouverEnd_A - thisMeanEndPos' );
    
end


%% Plots

figure(1)
hold on
title('Norm of Error between Rendezvous Point and Mean MC End Point')
plot( absDeviationEndPos )
hold off







