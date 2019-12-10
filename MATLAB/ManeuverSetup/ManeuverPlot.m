
run SatelliteInitialPositions2;

run earthParametersHPOP;


anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 2000;

numPeriods = 1;
numSamples = maneuverStartDelay + maneuverTime + 1;



%% Solve Lambert's problem to intersect new orbit

[rECIManeuverStart_target, vECIManeuverStart_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations )
[rECIManeuverStart_chaser, vECIManeuverStart_chaser] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations )


% Target's position after maneuver time
%[ rECIManeuverEnd_target, vECIManeuverEnd_target ] = nextStateTimeStep( muEarth, rECIManeuverStart_target, vECIManeuverStart_target, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations )
[ rECIManeuverEnd_target, vECIManeuverEnd_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverStartDelay + maneuverTime, anomalyErrorTolerance, anomalyMaxIterations )


% Required velocity change satellite B
[ deltaVStart_chaser, deltaVEnd_chaser, vIntersectOrbit ] = interceptOrbit( rECIManeuverStart_chaser, vECIManeuverStart_chaser, rECIManeuverEnd_target, vECIManeuverEnd_target, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations )

[ rECIManeuverEnd_chaser, vECIManeuverEnd_chaser ] = nextStateTimeStep( muEarth, rECIManeuverStart_chaser, vECIManeuverStart_chaser + deltaVStart_chaser, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations )

endPosError = rECIManeuverEnd_target - rECIManeuverEnd_chaser
endVelError = vIntersectOrbit - vECIManeuverEnd_chaser

%%


% [rECItest0_target1, vECItest0_target1] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, period_target*0.5, anomalyErrorTolerance, anomalyMaxIterations );
% [rECItest0_target2, vECItest0_target2] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, period_target*1.5, anomalyErrorTolerance, anomalyMaxIterations );
% endPosErrorTarget = rECItest0_target1 - rECItest0_target2
% 
% %[rECItest0_target, vECItest0_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, period_target, anomalyErrorTolerance, anomalyMaxIterations );
% %endPosErrorTarget = r0ECI_target - rECItest0_target
% 
% [rECItest1_target, vECItest1_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, period_target*0.75, anomalyErrorTolerance, anomalyMaxIterations );
% 
% [ deltaVStart_test2_target, deltaVEnd_test2_target, vIntersectOrbit_test2_target ] = interceptOrbit( r0ECI_target, v0ECI_target, rECItest1_target, vECItest1_target, period_target*0.5, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
% [rECItest2_target, vECItest2_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target + deltaVStart_test2_target, period_target*0.5, anomalyErrorTolerance, anomalyMaxIterations );
% endPosErrorTarget = rECItest1_target - rECItest2_target

%%

% test3_endPos_desired = [0; 7500; 0]
% [ deltaVStart_test3_target, deltaVEnd_test3_target, vIntersectOrbit_test3_target ] = interceptOrbit( r0ECI_target, v0ECI_target, test3_endPos_desired, 0, 3000, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations )
% [rECItest3_target, vECItest3_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target + deltaVStart_test3_target, 3000, anomalyErrorTolerance, anomalyMaxIterations )
% endPosErrorTarget = test3_endPos_desired - rECItest3_target



%% ECI motion

[rECI_targetX, rECI_targetY, rECI_targetZ, vECI_targetX, vECI_targetY, vECI_targetZ, sampleT_ECItarget] = ECITrajectory( r0ECI_target, v0ECI_target, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriods, numSamples, muEarth );

[rECI_chaserXinitial, rECI_chaserYinitial, rECI_chaserZinitial, vECI_chaserXinitial, vECI_chaserYinitial, vECI_chaserZinitial, sampleTinitial_ECIchaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriods, numSamples, muEarth );

if maneuverStartDelay == 0
    rECI_chaserX0 = r0ECI_chaser(1);
    rECI_chaserY0 = r0ECI_chaser(2);
    rECI_chaserZ0 = r0ECI_chaser(3);
    vECI_chaserX0 = v0ECI_chaser(1);
    vECI_chaserY0 = v0ECI_chaser(2);
    vECI_chaserZ0 = v0ECI_chaser(3);
else
    [rECI_chaserX0, rECI_chaserY0, rECI_chaserZ0, vECI_chaserX0, vECI_chaserY0, vECI_chaserZ0, sampleT0_ECIchaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
end
rECIcurr_chaser = rECIManeuverStart_chaser; %[ rECI_chaserX0( end ); rECI_chaserY0( end ); rECI_chaserZ0( end ) ];
vECIcurr_chaser = vECIManeuverStart_chaser; %[ vECI_chaserX0( end ); vECI_chaserY0( end ); vECI_chaserZ0( end ) ];
vECIcurr_chaser = vECIcurr_chaser;
[rECI_chaserX, rECI_chaserY, rECI_chaserZ, vECI_chaserX, vECI_chaserY, vECI_chaserZ, sampleT_ECIchaser] = ECITrajectory( rECIcurr_chaser, vECIcurr_chaser + deltaVStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );

[rECI_chaserX_unactuated, rECI_chaserY_unactuated, rECI_chaserZ_unactuated, vECI_chaserX_unactuated, vECI_chaserY_unactuated, vECI_chaserZ_unactuated, sampleT_ECIchaser_unactuated] = ECITrajectory( rECIcurr_chaser, vECIcurr_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );


%% Relative motion

if maneuverStartDelay == 0
    rLVLH_Rel = BPosRelativeToA( r0ECI_target, v0ECI_target, r0ECI_chaser );
    rLVLH_chaserX0 = rLVLH_Rel(1);
    rLVLH_chaserY0 = rLVLH_Rel(2);
    rLVLH_chaserZ0 = rLVLH_Rel(3);
else
    [rLVLH_chaserX0, rLVLH_chaserY0, rLVLH_chaserZ0, rLVLH_chaserNorm0, sampleTchaser0, lastECIPos_chaser0, lastECIVel_chaser0 ] = relativeTrajectory( r0ECI_target, v0ECI_target, r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
end
[rLVLH_chaserX1, rLVLH_chaserY1, rLVLH_chaserZ1, rLVLH_chaserNorm1, sampleTchaser1, lastECIPos_chaser1, lastECIVel_chaser1 ] = relativeTrajectory( rECIManeuverStart_target, vECIManeuverStart_target, rECIManeuverStart_chaser, vECIManeuverStart_chaser + deltaVStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );

[rLVLH_chaserX1_unactuated, rLVLH_chaserY1_unactuated, rLVLH_chaserZ1_unactuated, rLVLH_chaserNorm1_unactuated, sampleTchaser1_unactuated, lastECIPos_chaser1_unactuated, lastECIVel_chaser1_unactuated ] = relativeTrajectory( rECIManeuverStart_target, vECIManeuverStart_target, rECIManeuverStart_chaser, vECIManeuverStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );


%%

%rLVLH_chaserNorm1(end)

%%

if maneuverStartDelay == 0
    maneuverStartDelay = 1;
end


figure(1)
hold on
title('ECI Trajectories')
%axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
[ sx, sy, sz ] = sphere;
surf( sx*rEarth, sy*rEarth, sz*rEarth, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.05 );
h1 = plot3( rECI_targetX, rECI_targetY, rECI_targetZ, '-.');
h3 = plot3( rECI_chaserXinitial, rECI_chaserYinitial, rECI_chaserZinitial, '-.');
h2 = plot3( rECI_chaserX0, rECI_chaserY0, rECI_chaserZ0);
h4 = plot3( rECI_chaserX, rECI_chaserY, rECI_chaserZ);
% Draw the initial position vectors:
%plot3( [0 rECI_targetX( 1 ) ], [0 rECI_targetY( 1 ) ], [0 rECI_targetZ( 1 ) ], '-.' )
%plot3( [0 rECI_chaserX0(1) ], [0 rECI_chaserY0(1) ], [0 rECI_chaserZ0(1) ], '-.' )
%plot3( [0 rECI_chaserX( end ) ], [0 rECI_chaserY( end ) ], [0 rECI_chaserZ( end ) ], '-.' )
%plot3( [0 rECI_targetX( end ) ], [0 rECI_targetY( end ) ], [0 rECI_targetZ( end ) ], '-.' )
text( rECI_targetX( 1 ), rECI_targetY( 1 ), rECI_targetZ( 1 ), 'Target Initial Point' )
text( rECI_chaserX0(1), rECI_chaserY0(1), rECI_chaserZ0(1), 'Chaser Initial Point')
text( rECI_chaserX0( end ), rECI_chaserY0( end ), rECI_chaserZ0( end ), 'Chaser Maneuver Start' )
text( rECI_targetX( maneuverStartDelay ), rECI_targetY( maneuverStartDelay ), rECI_targetZ( maneuverStartDelay ), 'Target Maneuver Start' )
text( rECI_chaserX( end ), rECI_chaserY( end ), rECI_chaserZ( end ), 'Chaser Point of Rendezvous' )
text( rECI_targetX( end ), rECI_targetY( end ), rECI_targetZ( end ), 'Target Final Point' )
grid on
%legend([h1, h2, h3], 'Target Trajectory', 'Chaser Initial Trajectory', 'Chaser Maneuver Trajectory')
hold off


figure(2)
hold on
title('Experiment 1, 2 \& 3 Chaser Motion in Target LVLH')
plot3( rLVLH_chaserX0, rLVLH_chaserY0, rLVLH_chaserZ0, 'b' )
plot3( rLVLH_chaserX1, rLVLH_chaserY1, rLVLH_chaserZ1, 'g' )
plot3( rLVLH_chaserX1_unactuated, rLVLH_chaserY1_unactuated, rLVLH_chaserZ1_unactuated, 'b-.' )
plot3( 0, 0, 0, 'r*' )
axis equal
axis on
grid on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
% Label the origin of the moving frame attached to A:
%text (0, 0, 0, 'Target')
% Label the start of relative trajectories:
text(rLVLH_chaserX0(1), rLVLH_chaserY0(1), rLVLH_chaserZ0(1), 'Chaser Initial Position')
text(rLVLH_chaserX1( end ), rLVLH_chaserY1( end ), rLVLH_chaserZ1( end ), 'Chaser Final Position')
text(rLVLH_chaserX1(1), rLVLH_chaserY1(1), rLVLH_chaserZ1(1), 'Chaser Maneuver Start')
hold off




%%


% chaserManeuverStartPos = [rECI_chaserX0( end ), rECI_chaserY0( end ), rECI_chaserZ0( end )];
% targetManeuverStartPos = [rECI_targetX( maneuverStartDelay ), rECI_targetY( maneuverStartDelay ), rECI_targetZ( maneuverStartDelay )];
% 
% diffManeuverStartPos = targetManeuverStartPos - chaserManeuverStartPos
% 
% diffSimStartPos = r0ECI_target - r0ECI_chaser
% 
% targetEndPos = [rECI_targetX( end ), rECI_targetY( end ), rECI_targetZ( end )];
% chaserEndPos = [rECI_chaserX( end ), rECI_chaserY( end ), rECI_chaserZ( end )];
% 
% diffEndPos = targetEndPos - chaserEndPos

diffECItargetEnd = rECIManeuverEnd_target - [rECI_targetX( end ), rECI_targetY( end ), rECI_targetZ( end )]'

rendezvousErrorECIEnd = [rECI_chaserX( end ), rECI_chaserY( end ), rECI_chaserZ( end )]' - [rECI_targetX( end ), rECI_targetY( end ), rECI_targetZ( end )]'

