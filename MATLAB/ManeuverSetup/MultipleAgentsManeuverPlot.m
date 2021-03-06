
run SatelliteInitialPositions3;

run earthParametersHPOP;


anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 2000;

numPeriods = 1;
numSamples = maneuverStartDelay + maneuverTime + 1;



%% Solve Lambert's problem to intersect new orbit

[rECIManeuverStart_target, vECIManeuverStart_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations );
[rECIManeuverStart_chaser1, vECIManeuverStart_chaser1] = nextStateTimeStep( muEarth, r0ECI_chaser1, v0ECI_chaser1, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations );
[rECIManeuverStart_chaser2, vECIManeuverStart_chaser2] = nextStateTimeStep( muEarth, r0ECI_chaser2, v0ECI_chaser2, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations );
[rECIManeuverStart_chaser3, vECIManeuverStart_chaser3] = nextStateTimeStep( muEarth, r0ECI_chaser3, v0ECI_chaser3, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations );
[rECIManeuverStart_chaser4, vECIManeuverStart_chaser4] = nextStateTimeStep( muEarth, r0ECI_chaser4, v0ECI_chaser4, maneuverStartDelay, anomalyErrorTolerance, anomalyMaxIterations );


% Target's position after maneuver time
%[ rECIManeuverEnd_target, vECIManeuverEnd_target ] = nextStateTimeStep( muEarth, rECIManeuverStart_target, vECIManeuverStart_target, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations )
[ rECIManeuverEnd_target, vECIManeuverEnd_target ] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverStartDelay + maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );

QmatECItoLVLH_targetEnd = ECIToLVLH( rECIManeuverEnd_target, vECIManeuverEnd_target );
QmatLVLHtoECI_targetEnd = QmatECItoLVLH_targetEnd';

chaser1_ECI_endPos = rECIManeuverEnd_target + QmatLVLHtoECI_targetEnd*[0;0;chaserTargetDistance];
chaser2_ECI_endPos = rECIManeuverEnd_target + QmatLVLHtoECI_targetEnd*[chaserTargetDistance;0;0];
chaser3_ECI_endPos = rECIManeuverEnd_target + QmatLVLHtoECI_targetEnd*[0;0;-chaserTargetDistance];
chaser4_ECI_endPos = rECIManeuverEnd_target + QmatLVLHtoECI_targetEnd*[-chaserTargetDistance;0;0];

% Required velocity change satellite B
[ deltaVStart_chaser1, deltaVEnd_chaser1, vIntersectOrbit_chaser1 ] = interceptOrbit( rECIManeuverStart_chaser1, vECIManeuverStart_chaser1, chaser1_ECI_endPos, vECIManeuverEnd_target, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
[ deltaVStart_chaser2, deltaVEnd_chaser2, vIntersectOrbit_chaser2 ] = interceptOrbit( rECIManeuverStart_chaser2, vECIManeuverStart_chaser2, chaser2_ECI_endPos, vECIManeuverEnd_target, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
[ deltaVStart_chaser3, deltaVEnd_chaser3, vIntersectOrbit_chaser3 ] = interceptOrbit( rECIManeuverStart_chaser3, vECIManeuverStart_chaser3, chaser3_ECI_endPos, vECIManeuverEnd_target, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
[ deltaVStart_chaser4, deltaVEnd_chaser4, vIntersectOrbit_chaser4 ] = interceptOrbit( rECIManeuverStart_chaser4, vECIManeuverStart_chaser4, chaser4_ECI_endPos, vECIManeuverEnd_target, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

[ rECIManeuverEnd_chaser1, vECIManeuverEnd_chaser1 ] = nextStateTimeStep( muEarth, rECIManeuverStart_chaser1, vECIManeuverStart_chaser1 + deltaVStart_chaser1, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );
[ rECIManeuverEnd_chaser2, vECIManeuverEnd_chaser2 ] = nextStateTimeStep( muEarth, rECIManeuverStart_chaser2, vECIManeuverStart_chaser2 + deltaVStart_chaser2, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );
[ rECIManeuverEnd_chaser3, vECIManeuverEnd_chaser3 ] = nextStateTimeStep( muEarth, rECIManeuverStart_chaser3, vECIManeuverStart_chaser3 + deltaVStart_chaser3, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );
[ rECIManeuverEnd_chaser4, vECIManeuverEnd_chaser4 ] = nextStateTimeStep( muEarth, rECIManeuverStart_chaser4, vECIManeuverStart_chaser4 + deltaVStart_chaser4, maneuverTime, anomalyErrorTolerance, anomalyMaxIterations );


endPosError1 = rECIManeuverEnd_target - rECIManeuverEnd_chaser1
endVelError1 = vIntersectOrbit_chaser1 - vECIManeuverEnd_chaser1
endPosError2 = rECIManeuverEnd_target - rECIManeuverEnd_chaser2
endVelError2 = vIntersectOrbit_chaser2 - vECIManeuverEnd_chaser2
endPosError3 = rECIManeuverEnd_target - rECIManeuverEnd_chaser3
endVelError3 = vIntersectOrbit_chaser3 - vECIManeuverEnd_chaser3
endPosError4 = rECIManeuverEnd_target - rECIManeuverEnd_chaser4
endVelError4 = vIntersectOrbit_chaser4 - vECIManeuverEnd_chaser4


%% ECI motion

[rECI_targetX, rECI_targetY, rECI_targetZ, vECI_targetX, vECI_targetY, vECI_targetZ, sampleT_ECItarget] = ECITrajectory( r0ECI_target, v0ECI_target, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriods, numSamples, muEarth );

[rECI_chaser1_Xinitial, rECI_chaser1_Yinitial, rECI_chaser1_Zinitial, vECI_chaser1_Xinitial, vECI_chaser1_Yinitial, vECI_chaser1_Zinitial, sampleTinitial_ECIchaser] = ECITrajectory( r0ECI_chaser1, v0ECI_chaser1, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriods, numSamples, muEarth );
[rECI_chaser1_X0, rECI_chaser1_Y0, rECI_chaser1_Z0, vECI_chaser1_X0, vECI_chaser1_Y0, vECI_chaser1_Z0, sampleT0_ECIchaser1] = ECITrajectory( r0ECI_chaser1, v0ECI_chaser1, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
rECIcurr_chaser1 = rECIManeuverStart_chaser1;
vECIcurr_chaser1 = vECIManeuverStart_chaser1; 
[rECI_chaser1_X, rECI_chaser1_Y, rECI_chaser1_Z, vECI_chaser1_X, vECI_chaser1_Y, vECI_chaser1_Z, sampleT_ECIchaser1] = ECITrajectory( rECIcurr_chaser1, vECIcurr_chaser1 + deltaVStart_chaser1, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );

[rECI_chaser2_Xinitial, rECI_chaser2_Yinitial, rECI_chaser2_Zinitial, vECI_chaser2_Xinitial, vECI_chaser2_Yinitial, vECI_chaser2_Zinitial, sampleTinitial_ECIchaser] = ECITrajectory( r0ECI_chaser2, v0ECI_chaser2, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriods, numSamples, muEarth );
[rECI_chaser2_X0, rECI_chaser2_Y0, rECI_chaser2_Z0, vECI_chaser2_X0, vECI_chaser2_Y0, vECI_chaser2_Z0, sampleT0_ECIchaser2] = ECITrajectory( r0ECI_chaser2, v0ECI_chaser2, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
rECIcurr_chaser2 = rECIManeuverStart_chaser2;
vECIcurr_chaser2 = vECIManeuverStart_chaser2; 
[rECI_chaser2_X, rECI_chaser2_Y, rECI_chaser2_Z, vECI_chaser2_X, vECI_chaser2_Y, vECI_chaser2_Z, sampleT_ECIchaser2] = ECITrajectory( rECIcurr_chaser2, vECIcurr_chaser2 + deltaVStart_chaser2, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );

[rECI_chaser3_Xinitial, rECI_chaser3_Yinitial, rECI_chaser3_Zinitial, vECI_chaser3_Xinitial, vECI_chaser3_Yinitial, vECI_chaser3_Zinitial, sampleTinitial_ECIchaser] = ECITrajectory( r0ECI_chaser3, v0ECI_chaser3, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriods, numSamples, muEarth );
[rECI_chaser3_X0, rECI_chaser3_Y0, rECI_chaser3_Z0, vECI_chaser3_X0, vECI_chaser3_Y0, vECI_chaser3_Z0, sampleT0_ECIchaser3] = ECITrajectory( r0ECI_chaser3, v0ECI_chaser3, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
rECIcurr_chaser3 = rECIManeuverStart_chaser3;
vECIcurr_chaser3 = vECIManeuverStart_chaser3; 
[rECI_chaser3_X, rECI_chaser3_Y, rECI_chaser3_Z, vECI_chaser3_X, vECI_chaser3_Y, vECI_chaser3_Z, sampleT_ECIchaser3] = ECITrajectory( rECIcurr_chaser3, vECIcurr_chaser3 + deltaVStart_chaser3, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );

[rECI_chaser4_Xinitial, rECI_chaser4_Yinitial, rECI_chaser4_Zinitial, vECI_chaser4_Xinitial, vECI_chaser4_Yinitial, vECI_chaser4_Zinitial, sampleTinitial_ECIchaser] = ECITrajectory( r0ECI_chaser4, v0ECI_chaser4, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay + maneuverTime, numPeriods, numSamples, muEarth );
[rECI_chaser4_X0, rECI_chaser4_Y0, rECI_chaser4_Z0, vECI_chaser4_X0, vECI_chaser4_Y0, vECI_chaser4_Z0, sampleT0_ECIchaser4] = ECITrajectory( r0ECI_chaser4, v0ECI_chaser4, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
rECIcurr_chaser4 = rECIManeuverStart_chaser4;
vECIcurr_chaser4 = vECIManeuverStart_chaser4; 
[rECI_chaser4_X, rECI_chaser4_Y, rECI_chaser4_Z, vECI_chaser4_X, vECI_chaser4_Y, vECI_chaser4_Z, sampleT_ECIchaser4] = ECITrajectory( rECIcurr_chaser4, vECIcurr_chaser4 + deltaVStart_chaser4, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );



[rECI_chaser1_X_unactuated, rECI_chaser1_Y_unactuated, rECI_chaser1_Z_unactuated, vECI_chaser1_X_unactuated, vECI_chaser1_Y_unactuated, vECI_chaser1_Z_unactuated, sampleT_ECIchaser_unactuated] = ECITrajectory( rECIcurr_chaser1, vECIcurr_chaser1, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );


%% Relative motion

[rLVLH_chaser1_X0, rLVLH_chaser1_Y0, rLVLH_chaser1_Z0, rLVLH_chaser1_Norm0, sampleTchaser1_0, lastECIPos_chaser1_0, lastECIVel_chaser1_0 ] = relativeTrajectory( r0ECI_target, v0ECI_target, r0ECI_chaser1, v0ECI_chaser1, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
[rLVLH_chaser1_X1, rLVLH_chaser1_Y1, rLVLH_chaser1_Z1, rLVLH_chaser1_Norm1, sampleTchaser1_1, lastECIPos_chaser1_1, lastECIVel_chaser1_1 ] = relativeTrajectory( rECIManeuverStart_target, vECIManeuverStart_target, rECIManeuverStart_chaser1, vECIManeuverStart_chaser1 + deltaVStart_chaser1, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );

[rLVLH_chaser2_X0, rLVLH_chaser2_Y0, rLVLH_chaser2_Z0, rLVLH_chaser2_Norm0, sampleTchaser2_0, lastECIPos_chaser2_0, lastECIVel_chaser2_0 ] = relativeTrajectory( r0ECI_target, v0ECI_target, r0ECI_chaser2, v0ECI_chaser2, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
[rLVLH_chaser2_X1, rLVLH_chaser2_Y1, rLVLH_chaser2_Z1, rLVLH_chaser2_Norm1, sampleTchaser2_1, lastECIPos_chaser2_1, lastECIVel_chaser2_1 ] = relativeTrajectory( rECIManeuverStart_target, vECIManeuverStart_target, rECIManeuverStart_chaser2, vECIManeuverStart_chaser2 + deltaVStart_chaser2, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );

[rLVLH_chaser3_X0, rLVLH_chaser3_Y0, rLVLH_chaser3_Z0, rLVLH_chaser3_Norm0, sampleTchaser3_0, lastECIPos_chaser3_0, lastECIVel_chaser3_0 ] = relativeTrajectory( r0ECI_target, v0ECI_target, r0ECI_chaser3, v0ECI_chaser3, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
[rLVLH_chaser3_X1, rLVLH_chaser3_Y1, rLVLH_chaser3_Z1, rLVLH_chaser3_Norm1, sampleTchaser3_1, lastECIPos_chaser3_1, lastECIVel_chaser3_1 ] = relativeTrajectory( rECIManeuverStart_target, vECIManeuverStart_target, rECIManeuverStart_chaser3, vECIManeuverStart_chaser3 + deltaVStart_chaser3, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );

[rLVLH_chaser4_X0, rLVLH_chaser4_Y0, rLVLH_chaser4_Z0, rLVLH_chaser4_Norm0, sampleTchaser4_0, lastECIPos_chaser4_0, lastECIVel_chaser4_0 ] = relativeTrajectory( r0ECI_target, v0ECI_target, r0ECI_chaser4, v0ECI_chaser4, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartDelay, numPeriods, maneuverStartDelay + 1, muEarth );
[rLVLH_chaser4_X1, rLVLH_chaser4_Y1, rLVLH_chaser4_Z1, rLVLH_chaser4_Norm1, sampleTchaser4_1, lastECIPos_chaser4_1, lastECIVel_chaser4_1 ] = relativeTrajectory( rECIManeuverStart_target, vECIManeuverStart_target, rECIManeuverStart_chaser4, vECIManeuverStart_chaser4 + deltaVStart_chaser4, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );


[rLVLH_chaserX1_unactuated, rLVLH_chaserY1_unactuated, rLVLH_chaserZ1_unactuated, rLVLH_chaserNorm1_unactuated, sampleTchaser1_unactuated, lastECIPos_chaser1_unactuated, lastECIVel_chaser1_unactuated ] = relativeTrajectory( rECIManeuverStart_target, vECIManeuverStart_target, rECIManeuverStart_chaser1, vECIManeuverStart_chaser1, anomalyErrorTolerance, anomalyMaxIterations, maneuverTime, numPeriods, maneuverTime + 1, muEarth );


%%


figure(2)
hold on
set(gca,'FontSize',24)
set(gcf,'renderer','Painters','Position', [10 10 1300 500])
title('Chaser Maneuvers in Target LVLH frame')
%h1 = plot3( rLVLH_chaser1_X0, rLVLH_chaser1_Y0, rLVLH_chaser1_Z0,'LineWidth',1.5,'Color','r' )
h2 = plot3( rLVLH_chaser1_X1, rLVLH_chaser1_Y1, rLVLH_chaser1_Z1,'LineWidth',0.75)
h2.Color(4) = 0.5;
%h3 = plot3( rLVLH_chaser2_X0, rLVLH_chaser2_Y0, rLVLH_chaser2_Z0,'LineWidth',1.5,'Color','b' )
h4 = plot3( rLVLH_chaser2_X1, rLVLH_chaser2_Y1, rLVLH_chaser2_Z1,'LineWidth',0.75)
h4.Color(4) = 0.5;
%h5 = plot3( rLVLH_chaser3_X0, rLVLH_chaser3_Y0, rLVLH_chaser3_Z0,'LineWidth',1.5,'Color','g' )
h6 = plot3( rLVLH_chaser3_X1, rLVLH_chaser3_Y1, rLVLH_chaser3_Z1,'LineWidth',0.75 )
h6.Color(4) = 0.5;
%h7 = plot3( rLVLH_chaser4_X0, rLVLH_chaser4_Y0, rLVLH_chaser4_Z0,'LineWidth',1.5,'Color','c' )
h8 = plot3( rLVLH_chaser4_X1, rLVLH_chaser4_Y1, rLVLH_chaser4_Z1,'LineWidth',0.75)
h8.Color(4) = 0.5;
h9 = plot3( rLVLH_chaserX1_unactuated, rLVLH_chaserY1_unactuated, rLVLH_chaserZ1_unactuated,'LineWidth',0.75)
h9.Color(4) = 0.5;
h0 = plot3( 0, 0, 0, 'r*' )
axis equal
axis on
grid on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
legend([h0,h2,h4,h6,h8,h9],'Target', 'Chaser 1', 'Chaser 2', 'Chaser 3', 'Chaser 4', 'Unactuated Trajectory','Location','eastoutside')
% Label the origin of the moving frame attached to A:
%text (0, 0, 0, 'Target')
% Label the start of relative trajectories:
%text(rLVLH_chaserX0(1), rLVLH_chaserY0(1), rLVLH_chaserZ0(1)+20, 'Chaser Initial Position','FontSize', 15)
%text(rLVLH_chaserX1( end ), rLVLH_chaserY1( end ), rLVLH_chaserZ1( end )-40, 'Chaser Final Position','FontSize', 15)
%text(rLVLH_chaserX1(1), rLVLH_chaserY1(1)+20, rLVLH_chaserZ1(1), 'Chaser Maneuver Start','FontSize', 15)
view([165 10])
 ylim([-3.3 0.2])
 xlim([-0.9 0.2])
 zlim([-0.1 0.1])
hold off



