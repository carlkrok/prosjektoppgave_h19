% Initial parameters for spacecraft A and B. All based on deg, km and s.

run earthParameters; 
run satelliteParameters;

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 1000;


orbitType = "prograde";
orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );
numPeriods = 1;
numSamples = 1000;


manouverTime = 1000; % Seconds


deltaVScaleErrorVec = [ .95 : 0.01 : 1.05 ];


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



manouverEndPositionError = [ size(deltaVScaleErrorVec) ];
deltaVStartError = [ size(deltaVScaleErrorVec) ];
indexCounter = 1;

% A's position after manouver time
[ rECIManouverEnd_A, vECIManouverEnd_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, manouverTime, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStart_B, deltaVEnd_B, vIntersectOrbit ] = interceptOrbit( r0ECI_B, v0ECI_B, rECIManouverEnd_A, vECIManouverEnd_A, manouverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );


for deltaVScaleError = deltaVScaleErrorVec % [ percent ]
   
    
    deltaTimeAfterManouver = orbitPeriod_A - manouverTime;
    
    deltaVStartErronous_B = deltaVScaleError * deltaVStart_B;

    [rLVLH_RelB1X, rLVLH_RelB1Y, rLVLH_RelB1Z, rLVLH_RelB1Norm, sampleTB1, lastECIPos_B1, lastECIVel_B1 ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_B, v0ECI_B + deltaVStartErronous_B, anomalyErrorTolerance, anomalyMaxIterations, manouverTime, numPeriods, numSamples, muEarth );
    
    manouverEndPositionError( indexCounter ) = rLVLH_RelB1Norm( numSamples );
    deltaVStartError( indexCounter ) = sign( deltaVScaleError - 1 ) * norm( deltaVStartErronous_B - deltaVStart_B );
    indexCounter = indexCounter + 1;
    
end
    


figure(7)
hold on
title('Rendezvous Position Error with Reduced Thrust')
plot( deltaVScaleErrorVec, manouverEndPositionError )
xlabel('deltaVStart Scaling Error')
ylabel('Maneuver End Distance Error [km]')
axis on
grid on
hold off

figure(8)
hold on
title('Rendezvous Position Error with Reduced Thrust')
plot( deltaVStartError, manouverEndPositionError )
xlabel('deltaVStart Error [m/s]')
ylabel('Maneuver End Distance Error [km]')
axis on
grid on
hold off



























































% TODOs
% implement both analytical and numerical strategies for the propagation of the motion
% have a look at solution of the Kepler problem, 2.6 of the last book sent by Leonard (see also pages around 85), seems that equation 2.18 is
