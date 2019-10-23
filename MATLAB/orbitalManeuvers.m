% Initial parameters for spacecraft A and B. All based on deg, km and s.

run earthParameters; 
run satelliteParameters;


anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 1000;


orbitType = "prograde";
orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );
numPeriods = 1; %100;
numSamples = 1000; %10000;

maneuverTime = 2500; % Seconds
numPeriodsManeuver = 1;
numSamplesManeuver = 1000;


% Error variables on initial position, initial velocity and change of
% velocity

velChangeXScaleEroor = 0.99;
velChangeYScaleEroor = 0.99;
velChangeZScaleEroor = 0.99;

maneuverTimeDelay = 1; % [ seconds ]




%% Calculation of relative position, velocity and acceleration of
% spacecraft B relative to A

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


r0ECI_C = r0ECI_B;
r0ECI_D = r0ECI_B;
r0ECI_E = r0ECI_B;



%% Solve Lambert's problem to intersect new orbit

% A's position after delay time
[ rECIDelayEnd_A, vECIDelayEnd_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, maneuverTimeDelay, anomalyErrorTolerance, anomalyMaxIterations );

% A's position after maneuver time
[ rECIManeuverEnd_A, vECIManeuverEnd_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );

% A's position after delay and maneuver time
[ rECIDelayAndManeuverEnd_A, vECIDelayAndManeuverEnd_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, maneuverTimeDelay + maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStart_B, deltaVEnd_B, vIntersectOrbit ] = interceptOrbit( r0ECI_B, v0ECI_B, rECIManeuverEnd_A, vECIManeuverEnd_A, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

% Updating B and C's velocities
deltaVStartErronous_B = [deltaVStart_B(1)*velChangeXScaleEroor, deltaVStart_B(2)*velChangeYScaleEroor, deltaVStart_B(3)*velChangeZScaleEroor]';

vECI_B = v0ECI_B + deltaVStart_B;
vECI_C = v0ECI_B;
vECI_D = v0ECI_B + deltaVStartErronous_B;


%% Plot of relative motion

deltaTimeAfterManeuver = orbitPeriod_A - maneuverTime;

[rLVLH_RelB1X, rLVLH_RelB1Y, rLVLH_RelB1Z, rLVLH_RelB1Norm, sampleTB1, lastECIPos_B1, lastECIVel_B1 ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_B, vECI_B, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriodsManeuver, numSamplesManeuver, muEarth );
firstECIVel_B2 = lastECIVel_B1 + deltaVEnd_B;
[rLVLH_RelB2X, rLVLH_RelB2Y, rLVLH_RelB2Z, rLVLH_RelB2Norm, sampleTB2, lastECIPos_B2, lastECIVel_B2 ] = relativeTrajectory( rECIManeuverEnd_A, vECIManeuverEnd_A, lastECIPos_B1, firstECIVel_B2, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManeuver, numPeriods - numPeriodsManeuver, numSamples - numSamplesManeuver, muEarth );

[rLVLH_RelCX, rLVLH_RelCY, rLVLH_RelCZ, rLVLH_RelCNorm, sampleTC, lastECIPos_C, lastECIVel_C ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_C, vECI_C, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_A, numPeriods, numSamples, muEarth );

[rLVLH_RelD1X, rLVLH_RelD1Y, rLVLH_RelD1Z, rLVLH_RelD1Norm, sampleTD1, lastECIPos_D1, lastECIVel_D1 ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_B, vECI_D, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriodsManeuver, numSamplesManeuver, muEarth );
firstECIVel_D2 = lastECIVel_D1 + deltaVEnd_B;
[rLVLH_RelD2X, rLVLH_RelD2Y, rLVLH_RelD2Z, rLVLH_RelD2Norm, sampleTD2, lastECIPos_D2, lastECIVel_D2 ] = relativeTrajectory( rECIManeuverEnd_A, vECIManeuverEnd_A, lastECIPos_D1, firstECIVel_D2, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManeuver, numPeriods, numSamples, muEarth );

[rLVLH_RelE1X, rLVLH_RelE1Y, rLVLH_RelE1Z, rLVLH_RelE1Norm, sampleTE1, lastECIPos_E1, lastECIVel_E1 ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_B, v0ECI_B, anomalyErrorTolerance, anomalyMaxIterations, maneuverTimeDelay, numPeriodsManeuver, numSamplesManeuver, muEarth );
firstECIVel_E2 = lastECIVel_E1 + deltaVStart_B;
[rLVLH_RelE2X, rLVLH_RelE2Y, rLVLH_RelE2Z, rLVLH_RelE2Norm, sampleTE2, lastECIPos_E2, lastECIVel_E2 ] = relativeTrajectory( rECIDelayEnd_A, vECIDelayEnd_A, lastECIPos_E1, firstECIVel_E2, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriodsManeuver, numSamplesManeuver,muEarth );
firstECIVel_E3 = lastECIVel_E2 + deltaVEnd_B;
[rLVLH_RelE3X, rLVLH_RelE3Y, rLVLH_RelE3Z, rLVLH_RelE3Norm, sampleTE3, lastECIPos_E3, lastECIVel_E3 ] = relativeTrajectory( rECIDelayAndManeuverEnd_A, vECIDelayAndManeuverEnd_A, lastECIPos_E2, firstECIVel_E3, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManeuver - maneuverTimeDelay, numPeriods - numPeriodsManeuver, numSamples - numSamplesManeuver, muEarth );



[rECI_AX, rECI_AY, rECI_AZ, vECI_AX, vECI_AY, vECI_AZ, sampleT_ECIA] = ECITrajectory( r0ECI_A, v0ECI_A, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_A, numPeriods, numSamples, muEarth );

[rECI_B1X, rECI_B1Y, rECI_B1Z, vECI_B1X, vECI_B1Y, vECI_B1Z, sampleT_ECIB1] = ECITrajectory( r0ECI_B, vECI_B, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriodsManeuver, numSamplesManeuver, muEarth );
lastECIPos_B1 = [ rECI_B1X( numSamplesManeuver ), rECI_B1Y( numSamplesManeuver ), rECI_B1Z( numSamplesManeuver ) ]';
firstECIVel_B2 = [ vECI_B1X( numSamplesManeuver ), vECI_B1Y( numSamplesManeuver ), vECI_B1Z( numSamplesManeuver ) ]' + deltaVEnd_B;
[rECI_B2X, rECI_B2Y, rECI_B2Z, vECI_B2X, vECI_B2Y, vECI_B2Z, sampleT_ECIB2] = ECITrajectory( lastECIPos_B1 , firstECIVel_B2, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManeuver, numPeriods - numPeriodsManeuver, numSamples - numSamplesManeuver, muEarth );

[rECI_CX, rECI_CY, rECI_CZ, vECI_CX, vECI_CY, vECI_CZ, sampleT_ECIC] = ECITrajectory( r0ECI_B, v0ECI_B, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_A, numPeriods, numSamples, muEarth );

[rECI_D1X, rECI_D1Y, rECI_D1Z, vECI_D1X, vECI_D1Y, vECI_D1Z, sampleT_ECID1] = ECITrajectory( r0ECI_B, vECI_D, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriodsManeuver, numSamplesManeuver, muEarth );
lastECIPos_D1 = [ rECI_D1X( numSamplesManeuver ), rECI_D1Y( numSamplesManeuver ), rECI_D1Z( numSamplesManeuver ) ]';
firstECIVel_D2 = [ vECI_D1X( numSamplesManeuver ), vECI_D1Y( numSamplesManeuver ), vECI_D1Z( numSamplesManeuver ) ]' + deltaVEnd_B;
[rECI_D2X, rECI_D2Y, rECI_D2Z, vECI_D2X, vECI_D2Y, vECI_D2Z, sampleT_ECID2] = ECITrajectory( lastECIPos_D1 , firstECIVel_D2, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManeuver, numPeriods - numPeriodsManeuver, numSamples - numSamplesManeuver, muEarth );

[rECI_E1X, rECI_E1Y, rECI_E1Z, vECI_E1X, vECI_E1Y, vECI_E1Z, sampleT_ECIE1] = ECITrajectory( r0ECI_B, v0ECI_B, anomalyErrorTolerance, anomalyMaxIterations, maneuverTimeDelay, numPeriodsManeuver, numSamplesManeuver, muEarth );
lastECIPos_E1 = [ rECI_E1X( numSamplesManeuver ), rECI_E1Y( numSamplesManeuver ), rECI_E1Z( numSamplesManeuver ) ]';
firstECIVel_E2 = [ vECI_E1X( numSamplesManeuver ), vECI_E1Y( numSamplesManeuver ), vECI_E1Z( numSamplesManeuver ) ]' + deltaVStart_B;
[rECI_E2X, rECI_E2Y, rECI_E2Z, vECI_E2X, vECI_E2Y, vECI_E2Z, sampleT_ECIE2] = ECITrajectory( lastECIPos_E1, firstECIVel_E2, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriodsManeuver, numSamplesManeuver, muEarth );
lastECIPos_E2 = [ rECI_E2X( numSamplesManeuver ), rECI_E2Y( numSamplesManeuver ), rECI_E2Z( numSamplesManeuver ) ]';
firstECIVel_E3 = [ vECI_E2X( numSamplesManeuver ), vECI_E2Y( numSamplesManeuver ), vECI_E2Z( numSamplesManeuver ) ]' + deltaVEnd_B;
[rECI_E3X, rECI_E3Y, rECI_E3Z, vECI_E3X, vECI_E3Y, vECI_E3Z, sampleT_ECIE3] = ECITrajectory( lastECIPos_E2 , firstECIVel_E3, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManeuver - maneuverTimeDelay, numPeriods - numPeriodsManeuver, numSamples - numSamplesManeuver, muEarth );


%% Parameter error analysis










%%

[timeDifferenceIndexA, indexManeuverTimeA] = min(abs(sampleT_ECIA-maneuverTime));

figure(4)
hold on
title('Motion Relative to Satellite A')
plot3( [ rLVLH_RelB1X rLVLH_RelB2X ], [ rLVLH_RelB1Y rLVLH_RelB2Y ], [ rLVLH_RelB1Z rLVLH_RelB2Z ])
plot3( [ rLVLH_RelD1X rLVLH_RelD2X ], [ rLVLH_RelD1Y rLVLH_RelD2Y ], [ rLVLH_RelD1Z rLVLH_RelD2Z ])
plot3( [ rLVLH_RelE1X rLVLH_RelE2X rLVLH_RelE3X ], [ rLVLH_RelE1Y rLVLH_RelE2Y rLVLH_RelE3Y ], [ rLVLH_RelE1Z rLVLH_RelE2Z rLVLH_RelE3Z ])
plot3( rLVLH_RelCX, rLVLH_RelCY, rLVLH_RelCZ, '-' )
legend( 'B', 'B with thrust error', 'B with time delay', 'C' )
axis equal
axis on
grid on
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
% Label the origin of the moving frame attached to A:
text (0, 0, 0, 'A')
% Label the start of relative trajectories:
text(rLVLH_RelB1X(1), rLVLH_RelB1Y(1), rLVLH_RelB1Z(1), 'B')
text(rLVLH_RelCX(1), rLVLH_RelCY(1), rLVLH_RelCZ(1), 'C')
% Draw the initial position vectors:
%line([0 rLVLH_Rel1X(1)], [0 rLVLH_Rel1Y(1)], [0 rLVLH_Rel1Z(1)])
%line([0 rLVLH_Rel2X(1)], [0 rLVLH_Rel2Y(1)], [0 rLVLH_Rel2Z(1)])
hold off

figure(5)
hold on
title('Distance to Satellite A')
plot( [ sampleTB1 ( sampleTB2 + maneuverTime ) ], [ rLVLH_RelB1Norm rLVLH_RelB2Norm ] )
plot( [ sampleTD1 ( sampleTD2 + maneuverTime ) ], [ rLVLH_RelD1Norm rLVLH_RelD2Norm ] )
plot( [ sampleTE1 ( sampleTE2 + maneuverTimeDelay ) ( sampleTE3 + maneuverTime + maneuverTimeDelay ) ], [ rLVLH_RelE1Norm rLVLH_RelE2Norm rLVLH_RelE3Norm ] )
plot( sampleTC, rLVLH_RelCNorm )
legend( 'B', 'B with thrust error', 'B with time delay', 'C' )
xlabel('Time [s]')
ylabel('Distance [km]')
hold off


figure(6)
hold on
title('Rendezvous Maneuver')
[ sx, sy, sz ] = sphere;
surf( sx*rEarth, sy*rEarth, sz*rEarth, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.05 );
% plot3( rECI_AX, rECI_AY, rECI_AZ )
% plot3( rECI_B1X, rECI_B1Y, rECI_B1Z, 'Color','[0 1 0]' )
% plot3( rECI_B2X, rECI_B2Y, rECI_B2Z, 'Color','[0 1 0]' )
% plot3( rECI_D1X, rECI_D1Y, rECI_D1Z, 'Color','[1 0 0]' )
% plot3( rECI_D2X, rECI_D2Y, rECI_D2Z, 'Color','[1 0 0]' )
plot3( rECI_E1X, rECI_E1Y, rECI_E1Z, 'Color','[1 0 0]' )
plot3( rECI_E2X, rECI_E2Y, rECI_E2Z, 'Color','[1 0 0]' )
plot3( rECI_E3X, rECI_E3Y, rECI_E3Z, 'Color','[1 0 0]' )
% plot3( rECI_CX, rECI_CY, rECI_CZ, '--', 'Color','[0 0 1]' )
axis equal
% Draw the initial position vectors:
plot3( [0 rECI_AX( 1 ) ], [0 rECI_AY( 1 ) ], [0 rECI_AZ( 1 ) ], '-.' )
plot3( [0 rECI_B1X(1) ], [0 rECI_B1Y(1) ], [0 rECI_B1Z(1) ], '-.' )
plot3( [0 rECI_AX( indexManeuverTimeA ) ], [0 rECI_AY( indexManeuverTimeA ) ], [0 rECI_AZ( indexManeuverTimeA ) ], '-.' )
text( rECI_AX( 1 ), rECI_AY( 1 ), rECI_AZ( 1 ), 'A0' )
text( rECI_B1X(1), rECI_B1Y(1), rECI_B1Z(1), 'B0')
text( rECI_AX( indexManeuverTimeA ), rECI_AY( indexManeuverTimeA ), rECI_AZ( indexManeuverTimeA ), 'A at maneuverTime' )
axis on
grid on
hold off
