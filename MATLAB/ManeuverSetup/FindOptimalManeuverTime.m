
run SatelliteInitialPositions;

run earthParametersHPOP;

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 1000;


minTime = 180;
maxTime = 6000;
deltaTime = 60;

numpoints = ( maxTime - minTime ) / deltaTime;

firstManouver = [ numpoints ];
secondManouver = [ numpoints ];
totalDeltaV = [ numpoints ];

r0ECI_target = [ 6846.778547302488, -155.4912902730209, 1197.216618136735 ]';
v0ECI_target = [ -1.313213744521753, -0.961664408215407, 7.404405793820075 ]';
period_target = 5790.6;

manouverTimes = minTime : deltaTime : maxTime;

for thisManouverTime = manouverTimes
    
    index = ( ( thisManouverTime - minTime + deltaTime ) / deltaTime );
    
    [rECIEnd_target, vECIEnd_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, maneuverStartDelay + thisManouverTime, anomalyErrorTolerance, anomalyMaxIterations )
    
    [ deltaVStart, deltaVEnd, vIntersectOrbit ] = interceptOrbit( r0ECI_chaser, v0ECI_chaser, rECIEnd_target, vECIEnd_target, maneuverStartDelay + thisManouverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
    
    firstManouver( index ) = norm( deltaVStart );
    secondManouver( index ) = norm( deltaVEnd );
    totalDeltaV( index ) = firstManouver( index ) + secondManouver( index );
    
end

optimalManeuverTimes = islocalmin(totalDeltaV);

figure(1)
hold on
title('Velocity Change for Different Maneuver Times')
plot(minTime : deltaTime : maxTime, totalDeltaV)
plot(manouverTimes( optimalManeuverTimes ), totalDeltaV( optimalManeuverTimes ), '*')
xlabel('Maneuver Lengths [s]')
ylabel('Velocity Change [km/s]')
hold off



