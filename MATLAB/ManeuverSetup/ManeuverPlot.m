
run SatelliteInitialPositions;

run earthParametersHPOP;


anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 1000;


orbitType = "retrograde";

numPeriods = 1;
numSamples = 3000;

maneuverTime = 4680; % Seconds
numPeriodsManeuver = 1;
numSamplesManeuver = 1000;



%% Solve Lambert's problem to intersect new orbit

[rECIManeuverStart_target, vECIManeuverSart_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations );
[rECIManeuverStart_chaser, vECIManeuverSart_chaser] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations );


% Target's position after maneuver time
[ rECIManeuverEnd_target, vECIManeuverEnd_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverStartDelay + maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStart_chaser, deltaVEnd_chaser, vIntersectOrbit ] = interceptOrbit( rECIManeuverStart_chaser, vECIManeuverSart_chaser, rECIManeuverEnd_target, vECIManeuverEnd_target, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );


%% ECI motion

[rECI_targetX, rECI_targetY, rECI_targetZ, vECI_targetX, vECI_targetY, vECI_targetZ, sampleT_ECItarget] = ECITrajectory( r0ECI_target, v0ECI_target, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriods, numSamples, muEarth );

[rECI_chaserXinitial, rECI_chaserYinitial, rECI_chaserZinitial, vECI_chaserXinitial, vECI_chaserYinitial, vECI_chaserZinitial, sampleTinitial_ECIchaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriodsManeuver, numSamplesManeuver, muEarth );

[rECI_chaserX0, rECI_chaserY0, rECI_chaserZ0, vECI_chaserX0, vECI_chaserY0, vECI_chaserZ0, sampleT0_ECIchaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriodsManeuver, numSamplesManeuver, muEarth );
rECIcurr_chaser = [ rECI_chaserX0( end ); rECI_chaserY0( end ); rECI_chaserZ0( end ) ];
vECIcurr_chaser = [ vECI_chaserX0( end ); vECI_chaserY0( end ); vECI_chaserZ0( end ) ];
[rECI_chaserX, rECI_chaserY, rECI_chaserZ, vECI_chaserX, vECI_chaserY, vECI_chaserZ, sampleT_ECIchaser] = ECITrajectory( rECIcurr_chaser, vECIcurr_chaser + deltaVStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriodsManeuver, numSamplesManeuver, muEarth );


%% Relative motion

[rLVLH_chaserX0, rLVLH_chaserY0, rLVLH_chaserZ0, rLVLH_chaserNorm0, sampleTchaser0, lastECIPos_chaser0, lastECIVel_chaser0 ] = relativeTrajectory( r0ECI_target, v0ECI_target, r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriodsManeuver, numSamplesManeuver, muEarth );
[rLVLH_chaserX1, rLVLH_chaserY1, rLVLH_chaserZ1, rLVLH_chaserNorm1, sampleTchaser1, lastECIPos_chaser1, lastECIVel_chaser1 ] = relativeTrajectory( rECIManeuverStart_target, vECIManeuverSart_target, lastECIPos_chaser0, lastECIVel_chaser0 + deltaVStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriodsManeuver, numSamplesManeuver, muEarth );


%%


figure(1)
hold on
title('ECI Trajectories')
%axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
[ sx, sy, sz ] = sphere;
%surf( sx*rEarth, sy*rEarth, sz*rEarth, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.05 );
h1 = plot3( rECI_targetX, rECI_targetY, rECI_targetZ, '-.')
h3 = plot3( rECI_chaserXinitial, rECI_chaserYinitial, rECI_chaserZinitial, '-.')
h2 = plot3( rECI_chaserX0, rECI_chaserY0, rECI_chaserZ0)
h3 = plot3( rECI_chaserX, rECI_chaserY, rECI_chaserZ)
% Draw the initial position vectors:
plot3( [0 rECI_targetX( 1 ) ], [0 rECI_targetY( 1 ) ], [0 rECI_targetZ( 1 ) ], '-.' )
plot3( [0 rECI_chaserX0(1) ], [0 rECI_chaserY0(1) ], [0 rECI_chaserZ0(1) ], '-.' )
plot3( [0 rECI_chaserX( end ) ], [0 rECI_chaserY( end ) ], [0 rECI_chaserZ( end ) ], '-.' )
plot3( [0 rECI_targetX( end ) ], [0 rECI_targetY( end ) ], [0 rECI_targetZ( end ) ], '-.' )
text( rECI_targetX( 1 ), rECI_targetY( 1 ), rECI_targetZ( 1 ), 'Target Initial Point' )
text( rECI_chaserX0(1), rECI_chaserY0(1), rECI_chaserZ0(1), 'Chaser Initial Point')
text( rECI_chaserX( end ), rECI_chaserY( end ), rECI_chaserZ( end ), 'Point of Rendezvous' )
grid on
%legend([h1, h2, h3], 'Target Trajectory', 'Chaser Initial Trajectory', 'Chaser Maneuver Trajectory')
hold off


figure(2)
hold on
title('Chaser Motion Relative to Target')
plot3( [ rLVLH_chaserX0 rLVLH_chaserX1 ], [ rLVLH_chaserY0 rLVLH_chaserY1 ], [ rLVLH_chaserZ0 rLVLH_chaserZ1 ])
axis equal
axis on
grid on
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
% Label the origin of the moving frame attached to A:
text (0, 0, 0, 'Target')
% Label the start of relative trajectories:
text(rLVLH_chaserX0(1), rLVLH_chaserY0(1), rLVLH_chaserZ0(1), 'Chaser Initial Position')
text(rLVLH_chaserX1( end ), rLVLH_chaserY1( end ), rLVLH_chaserZ1( end ), 'Chaser Final Position')
text(rLVLH_chaserX1(1), rLVLH_chaserY1(1), rLVLH_chaserZ1(1), 'Chaser Maneuver Start')
hold off







