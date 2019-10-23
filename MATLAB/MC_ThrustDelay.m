%% Load Parameters

run earthParameters; 
run satelliteParameters;

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 2000;

orbitType = "prograde";
orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );

maneuverTime = 2000;


%% Monte Carlo Experiment Setup

MCsampleNum = 100;

meanDeviationTimeSetup = 0;
maxDeviationTimeSetup = 5;

MCdeviationTimes = meanDeviationTimeSetup - maxDeviationTimeSetup + ( 2 * maxDeviationTimeSetup * rand( MCsampleNum, 1 ) );


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
[ rECIExperimentStart_A, vECIExperimentStart_A ] = nextStateTimeStep( muEarth, r0ECI_A, v0ECI_A, 2 * maxDeviationTimeSetup, anomalyErrorTolerance, anomalyMaxIterations );
    
% A's position after ideal manouver
[ rECIManouverEnd_A, vECIManouverEnd_A ] = nextStateTimeStep( muEarth, rECIExperimentStart_A, vECIExperimentStart_A, manouverTime, anomalyErrorTolerance, anomalyMaxIterations );
 
% B's position at experiment start
[ rECIExperimentStart_B, vECIExperimentStart_B ] = nextStateTimeStep( muEarth, r0ECI_B, v0ECI_B, 2 * maxDeviationTimeSetup, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStartECI_B, deltaVEndECI_B, vIntersectOrbit_B ] = interceptOrbit( rECIExperimentStart_B, vECIExperimentStart_B, rECIManouverEnd_A, vECIManouverEnd_A, maneuverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );


%% Monte Carlo Experimetns










