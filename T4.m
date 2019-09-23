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
i_B = 0.01;
O_B = 90;
w_B = 90;
T_B = 0;

hNorm_C = 52059;
e_C = 0.025724;
i_C = 0.01;
O_C = 90;
w_C = 90;
T_C = 0;


anomalyErrorTolerance = 10^(-8);
anomalyMaxIterations = 1000;


orbitType = "prograde";
orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );
numPeriods = 1;
numSamples = 1000;

manouverTime = 1000; % Seconds


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

r0PQW_C = positionVectorPQW( muEarth, hNorm_C, e_C, T_C );
v0PQW_C = velocityVectorPQW( muEarth, hNorm_C, e_C, T_C );
QmatPQWtoECI_C = transformPQWtoECI( i_C, O_C, w_C );
r0ECI_C = QmatPQWtoECI_C * r0PQW_C;
v0ECI_C = QmatPQWtoECI_C * v0PQW_C;



%% Solve Lambert's problem to intersect new orbit


% A's position after manouver time

[ rECIManouverEnd_A, vECIManouverEnd_A ] = nextState( muEarth, r0ECI_A, v0ECI_A, manouverTime, anomalyErrorTolerance, anomalyMaxIterations );

% Required velocity change satellite B and C
[ deltaVStart_B, deltaVEnd_B, vIntersectOrbit ] = interceptOrbit( r0ECI_B, v0ECI_B, rECIManouverEnd_A, vECIManouverEnd_A, manouverTime, orbitType, muEarth, anomalyErrorTolerance, anomalyMaxIterations );


% Updating B and C's velocities 
vECI_B = v0ECI_B + deltaVStart_B;

vECI_C = v0ECI_C;



%% Plot of relative motion

[rLVLH_RelB1X, rLVLH_RelB1Y, rLVLH_RelB1Z, rLVLH_RelB1Norm, sampleTB1, lastECIPos_B1, lastECIVel_B1 ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_B, vECI_B, anomalyErrorTolerance, anomalyMaxIterations, manouverTime, numPeriods, numSamples, muEarth );

firstECIVel_B2 = lastECIVel_B1 + deltaVEnd_B;

[rLVLH_RelB2X, rLVLH_RelB2Y, rLVLH_RelB2Z, rLVLH_RelB2Norm, sampleTB2, lastECIPos_B2, lastECIVel_B2 ] = relativeTrajectory( rECIManouverEnd_A, vECIManouverEnd_A, lastECIPos_B1, firstECIVel_B2, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_A - manouverTime, numPeriods, numSamples, muEarth );

[rLVLH_RelCX, rLVLH_RelCY, rLVLH_RelCZ, rLVLH_RelCNorm, sampleTC, lastECIPos_C, lastECIVel_C ] = relativeTrajectory( r0ECI_A, v0ECI_A, r0ECI_C, vECI_C, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_A, numPeriods, numSamples, muEarth );


figure(4)
hold on
plot3( [ rLVLH_RelB1X rLVLH_RelB2X ], [ rLVLH_RelB1Y rLVLH_RelB2Y ], [ rLVLH_RelB1Z rLVLH_RelB2Z ], '-' )
plot3( rLVLH_RelCX, rLVLH_RelCY, rLVLH_RelCZ, '-' )
legend( 'B', 'C' )
%axis equal
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
plot( sampleTC, rLVLH_RelCNorm )
legend( 'B', 'C' )
hold off

%%




function [ deltaVStart, deltaVEnd, vIntersectOrbit ] = interceptOrbit( posStart, vStart, posEnd, vEnd, deltaTime, orbitTypeDebris, mu, anomalyTolerance, nMax )


    posStartNorm = norm( posStart );
    vStartNorm = norm( vStart );
    
    vStartRadial = dot( posStart, vStart ) / posStartNorm;
    vStartRadialNorm = norm( vStartRadial );
    
    alpha = 2 / posStartNorm - vStartNorm^2 / mu;
    
    X = findAnomaly( mu, alpha, posStartNorm, vStartRadialNorm, deltaTime, anomalyTolerance, nMax );
    
    f = fLagrange( posStartNorm, alpha, X );
    g = gLagrange( mu, alpha, X, deltaTime );
    
    r = f * posStart + g * vStart;
    rNorm = norm( r );
    
    df = dfLagrange( posStartNorm, rNorm, mu, alpha, X );
    dg = dgLagrange( rNorm, alpha, X );
    
    % Calculating required velocity in startPosition to arrive at end orbit
    % in deltaT seconds
    vRequiredStart = (1 / g) * (posEnd - f * posStart);
    vIntersectOrbit = (1 / g) * ( dg * posEnd - posStart );

    % Change in velocity required
    deltaVStart = vRequiredStart - vStart;
    deltaVEnd = vEnd - vIntersectOrbit;


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


function f = fLagrange( r0Norm, alpha, X )

    f = 1 - (X^2 / r0Norm) * C( alpha * X^2 );

end

function df = dfLagrange( r0Norm, rNorm, mu, alpha, X )

    df = sqrt( mu ) / (r0Norm * rNorm) * (alpha * X^3 * S( alpha * X^2 ) - X) ;

end

function dg = dgLagrange( rNorm, alpha, X )

    dg = 1 - (X^2 / rNorm) * C( alpha * X^2 );

end

function g = gLagrange( mu, alpha, X, deltaT )

    g = deltaT - (1 / sqrt( mu )) * X^3 * S( alpha * X^2 );

end


function X = findAnomaly( mu, alpha, r0Norm, v0RadialNorm, deltaT, tolerance, nMax )

    X = sqrt( mu ) * abs( alpha ) * deltaT;
    
    df = dF( r0Norm, v0RadialNorm, mu, X, alpha );
    f = F( r0Norm, v0RadialNorm, mu, X, alpha, deltaT );
    
    ratio = f / df;
    
    n = 0;
    while( abs(ratio) > tolerance && n < nMax )
        
        n = n + 1;
        df = dF( r0Norm, v0RadialNorm, mu, X, alpha );
        f = F( r0Norm, v0RadialNorm, mu, X, alpha, deltaT );
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

    
function val = F( r0Norm, v0RadialNorm, mu, X, alpha, deltaT )

    z = alpha * X^2;
    val = ((r0Norm * v0RadialNorm) / sqrt( mu )) * X^2 * C( z ) + (1 - alpha * r0Norm ) * X^3 * S( z ) + r0Norm * X - sqrt( mu ) * deltaT;
    
end

function val = dF( r0Norm, v0RadialNorm, mu, X, alpha )

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





