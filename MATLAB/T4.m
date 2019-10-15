
clear all
close all


%% Initial parameters for spacecraft A and B. All based on deg, km and s.

DEGCONV = pi/180;

rEarth = 6378;
muEarth = 398600;


hNorm_A = 52059;
e_A = 0.025724; % Eccentricity
i_A = 0; % Inclination
O_A = 0; % Right ascension
w_A = 0; % Argument of Perigee
T_A = 0; % True anomaly

hNorm_B = 52059;
e_B = 0.025724;
i_B = 15;
O_B = 0;
w_B = 30;
T_B = 0;

% hNorm_C = 52059;
% e_C = 0.025724;
% i_C = 15;
% O_C = 0;
% w_C = 30;
% T_C = 0;

% hNorm_D = 52059;
% e_D = 0.025724;
% i_D = 15;
% O_D = 0;
% w_D = 30;
% T_D = 0;

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 1000;


orbitType = "prograde";
orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );
numPeriods = 100;
numSamples = 10000;

manouverTime = 1000; % Seconds
numPeriodsManouver = 1;
numSamplesManouver = 1000;


% Error variables on initial position, initial velocity and change of
% velocity

velChangeXScaleEroor = 0.99;
velChangeYScaleEroor = 0.99;
velChangeZScaleEroor = 0.99;

manouverTimeDelay = 5; % [ seconds ]




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

% r0PQW_C = positionVectorPQW( muEarth, hNorm_C, e_C, T_C );
% v0PQW_C = velocityVectorPQW( muEarth, hNorm_C, e_C, T_C );
% QmatPQWtoECI_C = transformPQWtoECI( i_C, O_C, w_C );
% r0ECI_C = QmatPQWtoECI_C * r0PQW_C;
% v0ECI_C = QmatPQWtoECI_C * v0PQW_C;

% r0PQW_D = positionVectorPQW( muEarth, hNorm_D, e_D, T_D );
% v0PQW_D = velocityVectorPQW( muEarth, hNorm_D, e_D, T_D );
% QmatPQWtoECI_D = transformPQWtoECI( i_D, O_D, w_D );
% r0ECI_D = QmatPQWtoECI_D * r0PQW_D;
% v0ECI_D = QmatPQWtoECI_D * v0PQW_D;

r0ECI_C = r0ECI_B;
r0ECI_D = r0ECI_B;
r0ECI_E = r0ECI_B;



%% Solve Lambert's problem to intersect new orbit

% A's position after delay time
[ rECIDelayEnd_A, vECIDelayEnd_A ] = nextState( muEarth, r0ECI_A, v0ECI_A, manouverTimeDelay, anomalyErrorTolerance, anomalyMaxIterations );

% A's position after manouver time
[ rECIManouverEnd_A, vECIManouverEnd_A ] = nextState( muEarth, r0ECI_A, v0ECI_A, manouverTime, anomalyErrorTolerance, anomalyMaxIterations );

% A's position after delay and manouver time
[ rECIDelayAndManouverEnd_A, vECIDelayAndManouverEnd_A ] = nextState( muEarth, r0ECI_A, v0ECI_A, manouverTimeDelay + manouverTime, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B
[ deltaVStart_B, deltaVEnd_B, vIntersectOrbit ] = interceptOrbit( r0ECI_B, v0ECI_B, rECIManouverEnd_A, vECIManouverEnd_A, manouverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

% Updating B and C's velocities
deltaVStartErronous_B = [deltaVStart_B(1)*velChangeXScaleEroor, deltaVStart_B(2)*velChangeYScaleEroor, deltaVStart_B(3)*velChangeZScaleEroor]';

vECI_B = v0ECI_B + deltaVStart_B;
vECI_C = v0ECI_B;
vECI_D = v0ECI_B + deltaVStartErronous_B;


%% Plot of relative motion

deltaTimeAfterManouver = orbitPeriod_A - manouverTime;

[rLVLH_RelB1X, rLVLH_RelB1Y, rLVLH_RelB1Z, rLVLH_RelB1Norm, sampleTB1, lastECIPos_B1, lastECIVel_B1 ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_B, vECI_B, anomalyErrorTolerance, anomalyMaxIterations, manouverTime, numPeriodsManouver, numSamplesManouver, muEarth );
firstECIVel_B2 = lastECIVel_B1 + deltaVEnd_B;
[rLVLH_RelB2X, rLVLH_RelB2Y, rLVLH_RelB2Z, rLVLH_RelB2Norm, sampleTB2, lastECIPos_B2, lastECIVel_B2 ] = relativeTrajectory( rECIManouverEnd_A, vECIManouverEnd_A, lastECIPos_B1, firstECIVel_B2, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManouver, numPeriods - numPeriodsManouver, numSamples - numSamplesManouver, muEarth );

[rLVLH_RelCX, rLVLH_RelCY, rLVLH_RelCZ, rLVLH_RelCNorm, sampleTC, lastECIPos_C, lastECIVel_C ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_C, vECI_C, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_A, numPeriods, numSamples, muEarth );

[rLVLH_RelD1X, rLVLH_RelD1Y, rLVLH_RelD1Z, rLVLH_RelD1Norm, sampleTD1, lastECIPos_D1, lastECIVel_D1 ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_B, vECI_D, anomalyErrorTolerance, anomalyMaxIterations, manouverTime, numPeriodsManouver, numSamplesManouver, muEarth );
firstECIVel_D2 = lastECIVel_D1 + deltaVEnd_B;
[rLVLH_RelD2X, rLVLH_RelD2Y, rLVLH_RelD2Z, rLVLH_RelD2Norm, sampleTD2, lastECIPos_D2, lastECIVel_D2 ] = relativeTrajectory( rECIManouverEnd_A, vECIManouverEnd_A, lastECIPos_D1, firstECIVel_D2, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManouver, numPeriods, numSamples, muEarth );

[rLVLH_RelE1X, rLVLH_RelE1Y, rLVLH_RelE1Z, rLVLH_RelE1Norm, sampleTE1, lastECIPos_E1, lastECIVel_E1 ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_B, v0ECI_B, anomalyErrorTolerance, anomalyMaxIterations, manouverTimeDelay, numPeriodsManouver, numSamplesManouver, muEarth );
firstECIVel_E2 = lastECIVel_E1 + deltaVStart_B;
[rLVLH_RelE2X, rLVLH_RelE2Y, rLVLH_RelE2Z, rLVLH_RelE2Norm, sampleTE2, lastECIPos_E2, lastECIVel_E2 ] = relativeTrajectory( rECIDelayEnd_A, vECIDelayEnd_A, lastECIPos_E1, firstECIVel_E2, anomalyErrorTolerance, anomalyMaxIterations, manouverTime, numPeriodsManouver, numSamplesManouver,muEarth );
firstECIVel_E3 = lastECIVel_E2 + deltaVEnd_B;
[rLVLH_RelE3X, rLVLH_RelE3Y, rLVLH_RelE3Z, rLVLH_RelE3Norm, sampleTE3, lastECIPos_E3, lastECIVel_E3 ] = relativeTrajectory( rECIDelayAndManouverEnd_A, vECIDelayAndManouverEnd_A, lastECIPos_E2, firstECIVel_E3, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManouver - manouverTimeDelay, numPeriods - numPeriodsManouver, numSamples - numSamplesManouver, muEarth );



[rECI_AX, rECI_AY, rECI_AZ, vECI_AX, vECI_AY, vECI_AZ, sampleT_ECIA] = ECITrajectory( r0ECI_A, v0ECI_A, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_A, numPeriods, numSamples, muEarth );

[rECI_B1X, rECI_B1Y, rECI_B1Z, vECI_B1X, vECI_B1Y, vECI_B1Z, sampleT_ECIB1] = ECITrajectory( r0ECI_B, vECI_B, anomalyErrorTolerance, anomalyMaxIterations, manouverTime, numPeriodsManouver, numSamplesManouver, muEarth );
lastECIPos_B1 = [ rECI_B1X( numSamplesManouver ), rECI_B1Y( numSamplesManouver ), rECI_B1Z( numSamplesManouver ) ]';
firstECIVel_B2 = [ vECI_B1X( numSamplesManouver ), vECI_B1Y( numSamplesManouver ), vECI_B1Z( numSamplesManouver ) ]' + deltaVEnd_B;
[rECI_B2X, rECI_B2Y, rECI_B2Z, vECI_B2X, vECI_B2Y, vECI_B2Z, sampleT_ECIB2] = ECITrajectory( lastECIPos_B1 , firstECIVel_B2, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManouver, numPeriods - numPeriodsManouver, numSamples - numSamplesManouver, muEarth );

[rECI_CX, rECI_CY, rECI_CZ, vECI_CX, vECI_CY, vECI_CZ, sampleT_ECIC] = ECITrajectory( r0ECI_B, v0ECI_B, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_A, numPeriods, numSamples, muEarth );

[rECI_D1X, rECI_D1Y, rECI_D1Z, vECI_D1X, vECI_D1Y, vECI_D1Z, sampleT_ECID1] = ECITrajectory( r0ECI_B, vECI_D, anomalyErrorTolerance, anomalyMaxIterations, manouverTime, numPeriodsManouver, numSamplesManouver, muEarth );
lastECIPos_D1 = [ rECI_D1X( numSamplesManouver ), rECI_D1Y( numSamplesManouver ), rECI_D1Z( numSamplesManouver ) ]';
firstECIVel_D2 = [ vECI_D1X( numSamplesManouver ), vECI_D1Y( numSamplesManouver ), vECI_D1Z( numSamplesManouver ) ]' + deltaVEnd_B;
[rECI_D2X, rECI_D2Y, rECI_D2Z, vECI_D2X, vECI_D2Y, vECI_D2Z, sampleT_ECID2] = ECITrajectory( lastECIPos_D1 , firstECIVel_D2, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManouver, numPeriods - numPeriodsManouver, numSamples - numSamplesManouver, muEarth );

[rECI_E1X, rECI_E1Y, rECI_E1Z, vECI_E1X, vECI_E1Y, vECI_E1Z, sampleT_ECIE1] = ECITrajectory( r0ECI_B, v0ECI_B, anomalyErrorTolerance, anomalyMaxIterations, manouverTimeDelay, numPeriodsManouver, numSamplesManouver, muEarth );
lastECIPos_E1 = [ rECI_E1X( numSamplesManouver ), rECI_E1Y( numSamplesManouver ), rECI_E1Z( numSamplesManouver ) ]';
firstECIVel_E2 = [ vECI_E1X( numSamplesManouver ), vECI_E1Y( numSamplesManouver ), vECI_E1Z( numSamplesManouver ) ]' + deltaVStart_B;
[rECI_E2X, rECI_E2Y, rECI_E2Z, vECI_E2X, vECI_E2Y, vECI_E2Z, sampleT_ECIE2] = ECITrajectory( lastECIPos_E1, firstECIVel_E2, anomalyErrorTolerance, anomalyMaxIterations, manouverTime, numPeriodsManouver, numSamplesManouver, muEarth );
lastECIPos_E2 = [ rECI_E2X( numSamplesManouver ), rECI_E2Y( numSamplesManouver ), rECI_E2Z( numSamplesManouver ) ]';
firstECIVel_E3 = [ vECI_E2X( numSamplesManouver ), vECI_E2Y( numSamplesManouver ), vECI_E2Z( numSamplesManouver ) ]' + deltaVEnd_B;
[rECI_E3X, rECI_E3Y, rECI_E3Z, vECI_E3X, vECI_E3Y, vECI_E3Z, sampleT_ECIE3] = ECITrajectory( lastECIPos_E2 , firstECIVel_E3, anomalyErrorTolerance, anomalyMaxIterations, deltaTimeAfterManouver - manouverTimeDelay, numPeriods - numPeriodsManouver, numSamples - numSamplesManouver, muEarth );


%% Parameter error analysis










%%
indexManouverTimeA = round( numSamples * manouverTime / orbitPeriod_A );

figure(4)
hold on
plot3( [ rLVLH_RelB1X rLVLH_RelB2X ], [ rLVLH_RelB1Y rLVLH_RelB2Y ], [ rLVLH_RelB1Z rLVLH_RelB2Z ])
plot3( [ rLVLH_RelD1X rLVLH_RelD2X ], [ rLVLH_RelD1Y rLVLH_RelD2Y ], [ rLVLH_RelD1Z rLVLH_RelD2Z ])
plot3( [ rLVLH_RelE1X rLVLH_RelE2X rLVLH_RelE3X ], [ rLVLH_RelE1Y rLVLH_RelE2Y rLVLH_RelE3Y ], [ rLVLH_RelE1Z rLVLH_RelE2Z rLVLH_RelE3Z ])
plot3( rLVLH_RelCX, rLVLH_RelCY, rLVLH_RelCZ, '-' )
legend( 'B', 'B with thrust error', 'B with time delay', 'C' )
axis equal
axis on
grid on
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
plot( [ sampleTB1 ( sampleTB2 + manouverTime ) ], [ rLVLH_RelB1Norm rLVLH_RelB2Norm ] )
plot( [ sampleTD1 ( sampleTD2 + manouverTime ) ], [ rLVLH_RelD1Norm rLVLH_RelD2Norm ] )
plot( [ sampleTE1 ( sampleTE2 + manouverTimeDelay ) ( sampleTE3 + manouverTime + manouverTimeDelay ) ], [ rLVLH_RelE1Norm rLVLH_RelE2Norm rLVLH_RelE3Norm ] )
plot( sampleTC, rLVLH_RelCNorm )
legend( 'B', 'B with thrust error', 'B with time delay', 'C' )
hold off


figure(6)
hold on
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
plot3( [0 rECI_AX( indexManouverTimeA ) ], [0 rECI_AY( indexManouverTimeA ) ], [0 rECI_AZ( indexManouverTimeA ) ], '-.' )
text( rECI_AX( 1 ), rECI_AY( 1 ), rECI_AZ( 1 ), 'A0' )
text( rECI_B1X(1), rECI_B1Y(1), rECI_B1Z(1), 'B0')
text( rECI_AX( indexManouverTimeA ), rECI_AY( indexManouverTimeA ), rECI_AZ( indexManouverTimeA ), 'A at manouverTime' )
axis on
grid on
hold off



function [ deltaVStart, deltaVEnd, vIntersectOrbit ] = interceptOrbit( posStart, vStart, posEnd, vEnd, deltaTime, orbitTypeDebris, mu, tolerance, nMax )


   rStart = norm( posStart );
   rEnd = norm( posEnd );

   deltaTheta = acosd( dot( posStart, posEnd ) / (rStart * rEnd) );

   wOrbit = cross( posStart, posEnd );
   if ((wOrbit(3) < 0 && orbitTypeDebris == "prograde") || (~(wOrbit(3) < 0 ) && orbitTypeDebris == "retrograde"))
       deltaTheta = 360 - deltaTheta;
   end

   A = sind( deltaTheta ) * sqrt( (rStart * rEnd) / (1 - cosd( deltaTheta ) ) );

   % TODO: Find first estimate expression
   z0 = -100;
   while F( z0, rStart, rEnd, A, mu, deltaTime ) < 0
       z0 = z0 + 0.1;
   end
   z1 = newtonStep( z0, rStart, rEnd, A, mu, deltaTime );

   newtonCount = 0;
   while (abs( z1 - z0 ) > tolerance) && (newtonCount < nMax)

       z0 = z1;
       z1 = newtonStep( z0, rStart, rEnd, A, mu, deltaTime );
       newtonCount = newtonCount + 1;

   end

   z = z1;


   % Calculating Lagrange constants
   f = 1 - (y( z, rStart, rEnd, A ) / rStart);
   g = A * sqrt( y( z, rStart, rEnd, A ) / mu);
   df = (sqrt( mu ) / (rStart * rEnd)) * sqrt( y( z, rStart, rEnd, A ) / C( z ) ) * (z * S( z ) - 1);
   dg = 1 - (y( z, rStart, rEnd, A ) / rEnd);

   % Calculating required velocity in startPosition to arrive at end orbit
   % in deltaT seconds
   %...Equation 5.46a:
%...Equation 5.46b:
%...Equation 5.46d:
%...Equation 5.28:
%...Equation 5.29:
   vRequiredStart = (1 / g) * (posEnd - f * posStart);
   vIntersectOrbit = (1 / g) * ( dg * posEnd - posStart );

   % Change in velocity required
   deltaVStart = vRequiredStart - vStart;
   deltaVEnd = vEnd - vIntersectOrbit;

end


function val = y( z, r1, r2, A )

   val = r1 + r2 + A * ((z * S( z ) - 1) / sqrt( C( z ) ));
end


function z1 = newtonStep( z0, r1, r2, A, mu, deltaT )

   z1 = z0 - (F( z0, r1, r2, A, mu, deltaT ) / dF( z0, r1, r2, A ));

end


function val = dF( z, r1, r2, A )

   if z == 0
       val = (sqrt( 2 ) / 40) * y( 0, r1, r2, A )^(3 / 2) + (A / 8) * (sqrt( y( 0, r1, r2, A ) ) + A * sqrt( 1 / (2 * y( 0, r1, r2, A ) ) ));
   else
       val = (y( z, r1, r2, A ) / C( z ))^(3/2) * ((1 / (2 * z)) * (C( z ) - (3 / 2) * (S( z )^2 / C( z ))) + (3 / 4) * (S( z )^2 / C( z ))) + (A / 8) * (3 * (S(z ) / C( z )) * sqrt( y( z, r1, r2, A ) ) + A * sqrt( C( z ) / y( z, r1, r2, A )));
   end

end

function val = F( z, r1, r2, A, mu, deltaT )

   val = (y( z, r1, r2, A ) / C( z ))^(3 / 2) * S( z ) + A * sqrt( y( z, r1, r2, A ) ) - sqrt( mu ) * deltaT;

end

function [rLVLH_RelX, rLVLH_RelY, rLVLH_RelZ, rLVLH_RelNorm, sampleT, lastECIPos_B, lastECIVel_B ] = relativeTrajectory( r0_A, v0_A, r0_B, v0_B, anomalyTolerance, nMax, orbitPeriod_A, numPeriods, numSamples, muEarth )

   r_A = r0_A;
   v_A = v0_A;

   r_B = r0_B;
   v_B = v0_B;

   currTime = 0;
   endTime = currTime + orbitPeriod_A * numPeriods;
   sampleInterval = (orbitPeriod_A * numPeriods) / numSamples;

   r0LVLH_Rel = BPosRelativeToA( r0_A, v0_A, r0_B );

   rLVLH_RelX = [numSamples];
   rLVLH_RelY = [numSamples];
   rLVLH_RelZ = [numSamples];
   rLVLH_RelNorm = [numSamples];
   sampleT = [numSamples];

   rLVLH_RelX( 1 ) = r0LVLH_Rel( 1 );
   rLVLH_RelY( 1 ) = r0LVLH_Rel( 2 );
   rLVLH_RelZ( 1 ) = r0LVLH_Rel( 3 );
   rLVLH_RelNorm( 1 ) = norm(r0LVLH_Rel);
   sampleT( 1 ) = currTime;

   for sampleIter = 2 : 1 : numSamples

       currTime = currTime + sampleInterval;

       [r_A, v_A] = nextState( muEarth, r0_A, v0_A, currTime, anomalyTolerance, nMax );
       [r_B, v_B] = nextState( muEarth, r0_B, v0_B, currTime, anomalyTolerance, nMax );

       rLVLH_Rel = BPosRelativeToA( r_A, v_A, r_B );
       rLVLH_RelX( sampleIter ) = rLVLH_Rel( 1 );
       rLVLH_RelY( sampleIter ) = rLVLH_Rel( 2 );
       rLVLH_RelZ( sampleIter ) = rLVLH_Rel( 3 );
       rLVLH_RelNorm( sampleIter ) = norm(rLVLH_Rel);
       sampleT( sampleIter ) = currTime;

       if sampleIter == numSamples
           lastECIPos_B = r_B;
           lastECIVel_B = v_B;
       end

   end

end

function [rXECI, rYECI, rZECI, vXECI, vYECI, vZECI, sampleT] = ECITrajectory( r0, v0, anomalyTolerance, nMax, orbitPeriod_A, numPeriods, numSamples, muEarth )

   r = r0;
   v = v0;


   currTime = 0;
   endTime = currTime + orbitPeriod_A * numPeriods;
   sampleInterval = (orbitPeriod_A * numPeriods) / numSamples;

   rXECI = [numSamples];
   rYECI = [numSamples];
   rZECI = [numSamples];
   vXECI = [numSamples];
   vYECI = [numSamples];
   vZECI = [numSamples];
   sampleT = [numSamples];

   rXECI( 1 ) = r0( 1 );
   rYECI( 1 ) = r0( 2 );
   rZECI( 1 ) = r0( 3 );
   vXECI( 1 ) = v0( 1 );
   vYECI( 1 ) = v0( 2 );
   vZECI( 1 ) = v0( 3 );
   sampleT( 1 ) = currTime;

   for sampleIter = 2 : 1 : numSamples

       currTime = currTime + sampleInterval;

       [r, v] = nextState( muEarth, r0, v0, currTime, anomalyTolerance, nMax );

       rXECI( sampleIter ) = r( 1 );
       rYECI( sampleIter ) = r( 2 );
       rZECI( sampleIter ) = r( 3 );
       vXECI( sampleIter ) = v( 1 );
       vYECI( sampleIter ) = v( 2 );
       vZECI( sampleIter ) = v( 3 );
       sampleT( sampleIter ) = currTime;

   end

end

function [r, v] = nextState( mu, r0, v0, deltaT, anomalyTolerance, nMax )

   v0Radial = dot( r0, v0 ) / norm( r0 );
   v0RadialNorm = norm( v0Radial );

   r0Norm = norm( r0 );
   v0Norm = norm( v0 );

   alpha = 2 / r0Norm - v0Norm^2 / mu;

   X = findAnomaly( mu, alpha, r0Norm, v0RadialNorm, deltaT, anomalyTolerance, nMax );

   f = fLagrange( r0Norm, alpha, X );
   g = gLagrange( mu, alpha, X, deltaT );

   r = f * r0 + g * v0;
   rNorm = norm( r );

   df = dfLagrange( r0Norm, rNorm, mu, alpha, X );
   dg = dgLagrange( rNorm, alpha, X );

   v = df * r0 + dg * v0;

end


% Based on eq. 5.31 in [Curtis2011]
function f = fLagrange( r0Norm, alpha, X )

   f = 1 - (X^2 / r0Norm) * C( alpha * X^2 );

end

% Based on eq. 5.31 in [Curtis2011]
function df = dfLagrange( r0Norm, rNorm, mu, alpha, X )

   df = sqrt( mu ) / (r0Norm * rNorm) * (alpha * X^3 * S( alpha * X^2 ) - X) ;

end

% Based on eq. 5.31 in [Curtis2011]
function dg = dgLagrange( rNorm, alpha, X )

   dg = 1 - (X^2 / rNorm) * C( alpha * X^2 );

end

% Based on eq. 5.31 in [Curtis2011]
function g = gLagrange( mu, alpha, X, deltaT )

   g = deltaT - (1 / sqrt( mu )) * X^3 * S( alpha * X^2 );

end


function X = findAnomaly( mu, alpha, r0Norm, v0RadialNorm, deltaT, tolerance, nMax )

   X = sqrt( mu ) * abs( alpha ) * deltaT;

   df = anomalydF( r0Norm, v0RadialNorm, mu, X, alpha );
   f = anomalyF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT );

   ratio = f / df;

   n = 0;
   while( abs(ratio) > tolerance && n < nMax )

       n = n + 1;
       df = anomalydF( r0Norm, v0RadialNorm, mu, X, alpha );
       f = anomalyF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT );
       ratio = f / df;
       X = X - ratio;

   end

   if n == nMax
       fprintf('\n **No. iterations of Kepler''s equation = %g', n)
       fprintf('\n   f/df                                = %g\n', f/df)
   end

end


function rLVLH_Rel = BPosRelativeToA( rECI_A, vECI_A, rECI_B )

   hECI_A = cross( rECI_A, vECI_A );
   hECINorm_A = norm( hECI_A );

   i_unitVectorMovingFrame = ( rECI_A / norm( rECI_A ));
   k_unitVectorMovingFrame = ( hECI_A / hECINorm_A );
   j_unitVectorMovingFrame = cross( k_unitVectorMovingFrame, i_unitVectorMovingFrame );


   QmatECItoLVLH = [i_unitVectorMovingFrame';
                    j_unitVectorMovingFrame';
                    k_unitVectorMovingFrame'];

   rECI_Rel = rECI_B - rECI_A;

   rLVLH_Rel = QmatECItoLVLH * rECI_Rel;

end


function [rLVLH_Rel, vLVLH_Rel, aLVLH_Rel] = BRelativeToA( rECI_A, vECI_A, rECI_B, vECI_B )

   hECI_A = cross( rECI_A, vECI_A );
   hECINorm_A = norm( hECI_A );

   i_unitVectorMovingFrame = ( rECI_A / norm( rECI_A ));
   k_unitVectorMovingFrame = ( hECI_A / hECINorm_A );
   j_unitVectorMovingFrame = cross( k_unitVectorMovingFrame, i_unitVectorMovingFrame );


   QmatECItoLVLH = [i_unitVectorMovingFrame';
                    j_unitVectorMovingFrame';
                    k_unitVectorMovingFrame'];

   angularVelocityLVLH = hECI_A / (norm( rECI_A ))^2;
   angularAccelerationLVLH = -2 * (dot( vECI_A, rECI_A ) / (norm( rECI_A ))^2) * angularVelocityLVLH;

   aECI_A = -(muEarth / norm(rECI_A)^3) * rECI_A;
   aECI_B = -(muEarth / norm(rECI_B)^3) * rECI_B;


   rECI_Rel = rECI_B - rECI_A;
   vECI_Rel = vECI_B - vECI_A - cross( angularVelocityLVLH, rECI_Rel );
   aECI_Rel = aECI_B - aECI_A - cross( angularAccelerationLVLH, rECI_Rel ) - cross( angularVelocityLVLH, cross( angularVelocityLVLH, rECI_Rel )) - 2 * cross( angularVelocityLVLH, vECI_Rel );


   rLVLH_Rel = QmatECItoLVLH * rECI_Rel;
   vLVLH_Rel = QmatECItoLVLH * vECI_Rel;
   aLVLH_Rel = QmatECItoLVLH * aECI_Rel;

end


% Based on eq. 3.59 in [Curtis2011]
function val = anomalyF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT )

   z = alpha * X^2;
   val = ((r0Norm * v0RadialNorm) / sqrt( mu )) * X^2 * C( z ) + (1 - alpha * r0Norm ) * X^3 * S( z ) + r0Norm * X - sqrt( mu ) * deltaT;

end

% Based on eq. 3.64 in [Curtis2011]
function val = anomalydF( r0Norm, v0RadialNorm, mu, X, alpha )

   z = alpha * X^2;
   val = ((r0Norm * v0RadialNorm) / sqrt( mu )) * X * (1 - alpha * X^2 * S( z )) + (1 - alpha * r0Norm) * X^2 * C( z ) + r0Norm;

end


function val = S( z )

   if z > 0
       val = (sqrt( z ) - sin( sqrt( z ) )) / (sqrt( z ))^3;
   elseif z < 0
       val = (sinh( sqrt( -z ) ) - sqrt( -z )) / (sqrt( -z ))^3;
   else
       val = 1 / 6;
   end

end


function val = C( z )

   if z > 0
       val = (1 - cos( sqrt( z ) )) / z;
   elseif z < 0
       val = (cosh( sqrt( -z ) ) - 1) / -z;
   else
       val = 1 / 2;
   end

end


function T = orbitPeriod( mu, h, e )

   T  = (( 2 * pi ) / mu^2 ) * (h / sqrt( 1 - e^2 ))^3;

end


function Q = transformPQWtoECI( i, O, w )

   Q = [-sind( O ) * cosd( i ) * sind( w ) + cosd( O ) * cosd( w ), -sind( O ) * cosd( i ) * cosd( w ) - cosd( O ) * sind( w ),  sind( O ) * sind( i );
       cosd( O ) * cosd( i ) * sind( w ) + sind( O ) * cosd( w ),   cosd( O ) * cosd( i ) * cosd( w ) - sind( O ) * sind( w ),   -cosd( O ) * sind( i );
       sind( i ) * sind( w ),                                    sind( i ) * cosd( w ),                                   cosd( i )];

end


function r = positionVectorPQW( mu, h, e, Theta )

   r = ( h^2 / mu ) * ( 1 / ( 1 + e * cosd( Theta ))) * [ cosd( Theta ); sind( Theta ); 0 ];

end


function v = velocityVectorPQW( mu, h, e, Theta )

   v = ( mu / h ) * [ -sind( Theta ); e + cosd( Theta ); 0 ];

end





% TODOs
% implement both analytical and numerical strategies for the propagation of the motion
% have a look at solution of the Kepler problem, 2.6 of the last book sent by Leonard (see also pages around 85), seems that equation 2.18 is
