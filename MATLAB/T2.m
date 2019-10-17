% Initial parameters for spacecraft A and B. All based on deg, km and s.

run earthParameters;
run satelliteParameters;

% From example 7.2 in [Curtis2011]
% hNorm_A = 52059;
% e_A = 0.025724;
% i_A = 60;
% O_A = 40;
% w_A = 30;
% T_A = 40;
% hNorm_B = 52362;
% e_B = 0.0072696;
% i_B = 50;
% O_B = 40;
% w_B = 120;
% T_B = 40;

% Calculation of relative position, velocity and acceleration of
% spacecraft B relative to A

rPQW_A = positionVectorPQW( muEarth, hNorm_A, e_A, T_A );
vPQW_A = velocityVectorPQW( muEarth, hNorm_A, e_A, T_A );
QmatPQWtoECI_A = transformPQWtoECI( i_A, O_A, w_A );
rECI_A = QmatPQWtoECI_A * rPQW_A;
vECI_A = QmatPQWtoECI_A * vPQW_A;

rPQW_B = positionVectorPQW( muEarth, hNorm_B, e_B, T_B );
vPQW_B = velocityVectorPQW( muEarth, hNorm_B, e_B, T_B );
QmatPQWtoECI_B = transformPQWtoECI( i_B, O_B, w_B );
rECI_B = QmatPQWtoECI_B * rPQW_B;
vECI_B = QmatPQWtoECI_B * vPQW_B;

rPQW_C = positionVectorPQW( muEarth, hNorm_C, e_C, T_C );
vPQW_C = velocityVectorPQW( muEarth, hNorm_C, e_C, T_C );
QmatPQWtoECI_C = transformPQWtoECI( i_C, O_C, w_C );
rECI_C = QmatPQWtoECI_C * rPQW_C;
vECI_C = QmatPQWtoECI_C * vPQW_C;


%% Plot of relative motion

anomalyTolerance = 10^(-8);
nMax = 1000;

orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );
numPeriods = 1;
numSamples = 5000;

[rLVLH_Rel1X, rLVLH_Rel1Y, rLVLH_Rel1Z, rLVLH_Rel1Norm, sampleT1] = relativeTrajectory( rECI_A, vECI_A, rECI_B, vECI_B, anomalyTolerance, nMax, orbitPeriod_A, numPeriods, numSamples, muEarth );

[rLVLH_Rel2X, rLVLH_Rel2Y, rLVLH_Rel2Z, rLVLH_Rel2Norm, sampleT2] = relativeTrajectory( rECI_A, vECI_A, rECI_C, vECI_C, anomalyTolerance, nMax, orbitPeriod_A, numPeriods, numSamples, muEarth );


figure(1)
plot3( rLVLH_Rel1X, rLVLH_Rel1Y, rLVLH_Rel1Z, '-', rLVLH_Rel2X, rLVLH_Rel2Y, rLVLH_Rel2Z, '-' )
hold on
%axis equal
axis on
grid on
% Label the origin of the moving frame attached to A:
text (0, 0, 0, 'A')
% Label the start of relative trajectories:
text(rLVLH_Rel1X(1), rLVLH_Rel1Y(1), rLVLH_Rel1Z(1), 'B')
text(rLVLH_Rel2X(1), rLVLH_Rel2Y(1), rLVLH_Rel2Z(1), 'C')
% Draw the initial position vectors:
line([0 rLVLH_Rel1X(1)], [0 rLVLH_Rel1Y(1)], [0 rLVLH_Rel1Z(1)])
line([0 rLVLH_Rel2X(1)], [0 rLVLH_Rel2Y(1)], [0 rLVLH_Rel2Z(1)])
hold off

figure(2)
hold on
plot( sampleT1, rLVLH_Rel1Norm)
plot( sampleT2, rLVLH_Rel2Norm)
legend('B', 'C')
hold off


% TODOs
% implement both analytical and numerical strategies for the propagation of the motion
% have a look at solution of the Kepler problem, 2.6 of the last book sent by Leonard (see also pages around 85), seems that equation 2.18 is 

