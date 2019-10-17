% Initial parameters for spacecraft A and B. All based on deg, km and s.

run earthParameters; 
run satelliteParameters;

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 1000;


orbitType = "prograde";
orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );
numPeriods = 1;
numSamples = 5000;

multipleOrbitsNumPeriods = 20;
multipleOrbitsNumSamples = 10000;
multipleOrbitsManouverTimeDelay = 1;

multipleOrbitsSampleTimeVec = 0 : ( multipleOrbitsNumPeriods * orbitPeriod_A ) / multipleOrbitsNumSamples : multipleOrbitsNumPeriods * orbitPeriod_A;

manouverTime = 1000; % Seconds

manouverTimeDelayVec = [ 0 : 0.1 : 10 ]; % Seconds



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

r0ECI_E = r0ECI_B;


multipleOrbitsPositionError = [ multipleOrbitsNumSamples ];

manouverEndPositionError = [ size(manouverTimeDelayVec) ];



% A's position after manouver time
[ rECIManouverEnd_A, vECIManouverEnd_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, manouverTime, anomalyErrorTolerance, anomalyMaxIterations );
    
% Required velocity change satellite B
[ deltaVStartECI_B, deltaVEndECI_B, vIntersectOrbit ] = interceptOrbit( r0ECI_B, v0ECI_B, rECIManouverEnd_A, vECIManouverEnd_A, manouverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

[ QmatECItoLVLH1 ] = ECIToLVLH( r0ECI_B, v0ECI_B );
deltaVStartLVLH_B = QmatECItoLVLH1 * deltaVStartECI_B;
deltaVEndLVLH_B = QmatECItoLVLH1 * deltaVEndECI_B;

indexCounter = 1;
[ rECIDelayAndManouverEnd_A, vECIDelayAndManouverEnd_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, multipleOrbitsManouverTimeDelay + manouverTime, anomalyErrorTolerance, anomalyMaxIterations );
[ rECIDelayEnd_B, vECIDelayEnd_B ] = nextStateTimeStep( muEarth, r0ECI_B, v0ECI_B, multipleOrbitsManouverTimeDelay, anomalyErrorTolerance, anomalyMaxIterations );
[ QmatECItoLVLH1 ] = ECIToLVLH( rECIDelayEnd_B, vECIDelayEnd_B );
QmatLVLHtoECI1 = QmatECItoLVLH1';
deltaVStart_B1 = QmatLVLHtoECI1 * deltaVStartLVLH_B;
[ rECIDelayAndManouverEnd_B, vECIDelayAndManouverEnd_B ] = nextStateTimeStep( muEarth, rECIDelayEnd_B, vECIDelayEnd_B + deltaVStartLVLH_B, manouverTime, anomalyErrorTolerance, anomalyMaxIterations );
[ QmatECItoLVLH2 ] = ECIToLVLH( rECIDelayAndManouverEnd_B, vECIDelayAndManouverEnd_B );
QmatLVLHtoECI2 = QmatECItoLVLH2';
deltaVEnd_B2 = QmatLVLHtoECI2 * deltaVEndLVLH_B;
for sampleTime = multipleOrbitsSampleTimeVec % [ seconds ]

    % A's position after delay and manouver time
    [ rECIFinal_A, vECIFinal_A ] = nextStateTimeStep( muEarth, rECIDelayAndManouverEnd_A, vECIDelayAndManouverEnd_A, sampleTime, anomalyErrorTolerance, anomalyMaxIterations );

    
    [ rECIFinal_B, vECIFinal_B ] = nextStateTimeStep( muEarth, rECIDelayAndManouverEnd_B, vECIDelayAndManouverEnd_B + deltaVEnd_B2, sampleTime, anomalyErrorTolerance, anomalyMaxIterations );
    
    multipleOrbitsPositionError( indexCounter ) = norm( rECIFinal_A - rECIFinal_B );
    indexCounter = indexCounter + 1;
    
end
    

indexCounter = 1;
for manouverTimeDelay = manouverTimeDelayVec % [ seconds ]
    
    % A's position after delay and manouver time
    [ rECIDelayAndManouverEnd_A, vECIDelayAndManouverEnd_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, manouverTimeDelay + manouverTime, anomalyErrorTolerance, anomalyMaxIterations );


    [ rECIDelayEnd_B, vECIDelayEnd_B ] = nextStateTimeStep( muEarth, r0ECI_B, v0ECI_B, manouverTimeDelay, anomalyErrorTolerance, anomalyMaxIterations );
    
    [ QmatECItoLVLH3 ] = ECIToLVLH( rECIDelayEnd_B, vECIDelayEnd_B );
    QmatLVLHtoECI3 = QmatECItoLVLH3';
    deltaVStart_B3 = QmatLVLHtoECI3 * deltaVStartLVLH_B;
    
    [ rECIDelayAndManouverEnd_B, vECIDelayAndManouverEnd_B ] = nextStateTimeStep( muEarth, rECIDelayEnd_B, vECIDelayEnd_B + deltaVStart_B3, manouverTimeDelay + manouverTime, anomalyErrorTolerance, anomalyMaxIterations );
    
    manouverEndPositionError( indexCounter ) = norm( rECIDelayAndManouverEnd_A - rECIDelayAndManouverEnd_B );
    indexCounter = indexCounter + 1;
    
end



figure(9)
hold on
plot( manouverEndPositionError )
axis equal
grid on
hold off

figure(10)
hold on
plot( multipleOrbitsPositionError )
grid on
hold off

