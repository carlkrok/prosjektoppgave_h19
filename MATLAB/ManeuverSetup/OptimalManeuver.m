
run SatelliteInitialPositions3;

run earthParametersHPOP;

TOLERANCE = 10^(-7);

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 3000;

MAXDELTAV = 10;

maneuverLengths = 60 : 60 : period_target;

maneuverStartTimes = 0 : 60 : period_target;

numLengths = length( maneuverLengths );
numStartTimes = length( maneuverStartTimes );

firstManouver = zeros( numStartTimes, numLengths );
secondManouver = zeros( numStartTimes, numLengths );
totalDeltaV = zeros( numStartTimes, numLengths );
optimalManeuverTimes = zeros( numStartTimes, numLengths );

for thisManeuverStartTimeIndex = 1 : numStartTimes
    
    thismaneuverStartTime = maneuverStartTimes( thisManeuverStartTimeIndex );
    
    [rStartECI_target, vStartECI_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, thismaneuverStartTime, anomalyErrorTolerance, anomalyMaxIterations );
    [rStartECI_chaser, vStartECI_chaser] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, thismaneuverStartTime, anomalyErrorTolerance, anomalyMaxIterations );
    
    for thisManouverLengthIndex = 1 : numLengths
        
        thisManouverLength = maneuverLengths( thisManouverLengthIndex );

        [rECIEnd_target, vECIEnd_target] = nextStateTimeStep( muEarth, rStartECI_target, vStartECI_target, thisManouverLength, anomalyErrorTolerance, anomalyMaxIterations );

        [ deltaVStart, deltaVEnd, vIntersectOrbit ] = interceptOrbit( rStartECI_chaser, vStartECI_chaser, rECIEnd_target, vECIEnd_target,  thisManouverLength, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

        [rECIEnd_chaser, vECIEnd_chaser] = nextStateTimeStep( muEarth, rStartECI_chaser, vStartECI_chaser + deltaVStart, thisManouverLength, anomalyErrorTolerance, anomalyMaxIterations );

            
        if ( norm( deltaVStart ) > MAXDELTAV ) 
            deltaVStart = [MAXDELTAV, 0, 0];
        end
        
        if ( norm( deltaVEnd ) > MAXDELTAV ) 
            deltaVEnd = [MAXDELTAV, 0, 0];
        end
        
        [rECIEndRef_target, vECIEndRef_target] = nextStateTimeStep( muEarth, r0ECI_target, v0ECI_target, thismaneuverStartTime + thisManouverLength, anomalyErrorTolerance, anomalyMaxIterations );

        if norm( rECIEndRef_target - rECIEnd_target ) > TOLERANCE
            deltaVStart = [MAXDELTAV, 0, 0];
            deltaVEnd = [MAXDELTAV, 0, 0];
        end
        
        if norm( rECIEndRef_target - rECIEnd_chaser ) > TOLERANCE
            deltaVStart = [MAXDELTAV, 0, 0];
            deltaVEnd = [MAXDELTAV, 0, 0];
        end
        
        if norm( rECIEnd_target - rECIEnd_chaser ) > TOLERANCE
            deltaVStart = [MAXDELTAV, 0, 0];
            deltaVEnd = [MAXDELTAV, 0, 0];
        end

        firstManouver( thisManeuverStartTimeIndex, thisManouverLengthIndex ) = norm( deltaVStart );
        secondManouver( thisManeuverStartTimeIndex, thisManouverLengthIndex ) = norm( deltaVEnd );
        totalDeltaV( thisManeuverStartTimeIndex, thisManouverLengthIndex ) = firstManouver( thisManeuverStartTimeIndex, thisManouverLengthIndex ) + secondManouver( thisManeuverStartTimeIndex, thisManouverLengthIndex );

    end
    
    optimalManeuverTimes( thisManeuverStartTimeIndex, : ) = islocalmin( totalDeltaV( thisManeuverStartTimeIndex, : ) );
    
end

%%

figure(1)
hold on
title('Velocity Change for Different Maneuver Start Times and Lengths')
for thisManeuverStartTimeIndex = 1 : numStartTimes
    plot(maneuverLengths, totalDeltaV( thisManeuverStartTimeIndex, : ) )
    for optimalIndexInter = 1 : numLengths
        isOptimalIndex = optimalManeuverTimes( thisManeuverStartTimeIndex, optimalIndexInter );
        if isOptimalIndex > 0
            plot(maneuverLengths( optimalIndexInter ), totalDeltaV( thisManeuverStartTimeIndex, optimalIndexInter), '*')
        end
    end
end
xlabel('Maneuver Lengths [s]')
ylabel('Velocity Change [km/s]')
hold off

%%

surfX = zeros(1, numLengths * numStartTimes);
surfY = zeros(1, numLengths * numStartTimes);
surfZ = zeros(1, numLengths * numStartTimes);


% X is Start times
% Y is Maneuver times
% Z is Total deltaV

surfY(1 , 1 : numLengths ) = maneuverLengths;
for thisManeuverLengthIndex = 1 : (numStartTimes - 1)
    surfY( 1 , ( thisManeuverLengthIndex * numLengths ) + 1 : ( ( thisManeuverLengthIndex + 1 ) * numLengths ) ) = maneuverLengths;
end

for thisManeuverStartTimeIndex = 1 : numStartTimes
    thisManeuverStartTime = maneuverStartTimes( thisManeuverStartTimeIndex );
    surfX( 1 , ( (thisManeuverStartTimeIndex - 1) * numLengths ) + 1 : ( ( thisManeuverStartTimeIndex ) * numLengths ) ) = ones( 1, numLengths ) .* thisManeuverStartTime;
end

for thisManeuverStartTimeIndex = 1 : numStartTimes
    thisManeuverStartTime = maneuverStartTimes( thisManeuverStartTimeIndex );
    for thisManeuverLengthIndex = 1 : (numStartTimes - 1)
        surfZ( 1, (( thisManeuverStartTimeIndex - 1 ) *  numLengths ) + thisManeuverLengthIndex ) = totalDeltaV( thisManeuverStartTimeIndex, thisManeuverLengthIndex );
    end
end

%%


minTotalVector = min(totalDeltaV);

minTotDeltaV = min(minTotalVector)
[minTotRow,minTotCol] = find(totalDeltaV==minTotDeltaV)

minStartVector = min(firstManouver); 

minStartDeltaV = min(minStartVector)
[minStartRow,minStartCol] = find(firstManouver==minStartDeltaV)

%%

% maneuverStartTimes(minTotCol)
% maneuverLengths(minTotRow)
maneuverStartTimes(minTotRow )
maneuverLengths( minTotCol)

maneuverStartTimes(minStartRow)
maneuverLengths(minStartCol)

figure(2)
hold on
title('Total Velocity Change for Different Maneuver Start Times and Lengths')
surf(maneuverLengths, maneuverStartTimes, totalDeltaV)
%plot3(maneuverStartTimes(minTotCol), maneuverLengths(minTotRow), minTotDeltaV, '*')
plot3(maneuverLengths(minTotCol), maneuverStartTimes(minTotRow), minTotDeltaV, '*')
xlabel('Maneuver Start Times [s]')
ylabel('Maneuver Lengths [s]')
zlabel('Velocity Change [km/s]')
hold off

figure(3)
hold on
title('Initial Velocity Change for Different Maneuver Start Times and Lengths')
surf(maneuverLengths, maneuverStartTimes, firstManouver)
plot3(maneuverLengths(minStartCol), maneuverStartTimes(minStartRow), minStartDeltaV, '*')
xlabel('Maneuver Start Times [s]')
ylabel('Maneuver Lengths [s]')
zlabel('Velocity Change [km/s]')
hold off

