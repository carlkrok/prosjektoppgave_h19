%% Initial parameters for spacecraft A and B. All based on deg, km and s.

rEarth = 6378; 
muEarth = 398600;

hNorm_A = 52059;
e_A = 0.025724;
i_A = 60;
O_A = 40;
w_A = 30;
T_A = 40;

hNorm_B = 52362;
e_B = 0.0072696;
i_B = 59.9;
O_B = 40;
w_B = 120;
T_B = 40;


%% Calculation of relative position, velocity and acceleration of
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


%% Plot of relative motion

r0_A = rECI_A;
v0_A = vECI_A;

r_A = r0_A;
v_A = v0_A;

r0_B = rECI_B;
v0_B = vECI_B;

r_B = r0_B;
v_B = v0_B;

currTime = 0;
numSamples = 5000;
anomalyTolerance = 10^(-6);
nMax = 1000;

orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A );
numPeriods = 60;
endTime = currTime + orbitPeriod_A * numPeriods;
sampleInterval = (orbitPeriod_A * numPeriods) / numSamples;

r0LVLH_Rel = BPosRelativeToA( r0_A, v0_A, r0_B );

rLVLH_RelX = [numSamples];
rLVLH_RelY = [numSamples];
rLVLH_RelZ = [numSamples];

rLVLH_RelX( 1 ) = r0LVLH_Rel( 1 );
rLVLH_RelY( 1 ) = r0LVLH_Rel( 2 );
rLVLH_RelZ( 1 ) = r0LVLH_Rel( 3 );

for sampleIter = 2 : 1 : numSamples
   
    [r_A, v_A] = nextState( muEarth, r0_A, v0_A, currTime, anomalyTolerance, nMax );
    [r_B, v_B] = nextState( muEarth, r0_B, v0_B, currTime, anomalyTolerance, nMax );
    
    rLVLH_Rel = BPosRelativeToA( r_A, v_A, r_B );
    rLVLH_RelX( sampleIter ) = rLVLH_Rel( 1 );
    rLVLH_RelY( sampleIter ) = rLVLH_Rel( 2 );
    rLVLH_RelZ( sampleIter ) = rLVLH_Rel( 3 );
    
    currTime = currTime + sampleInterval;

end

figure(1)
plot3( rLVLH_RelX, rLVLH_RelY, rLVLH_RelZ, '*' )


%%



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
    
    f = F( r0Norm, v0RadialNorm, mu, X, alpha );
    df = dF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT );
    
    ratio = f / df;
    
    n = 0;
    while( abs(ratio) > tolerance && n < nMax )
        
        n = n + 1;
        f = F( r0Norm, v0RadialNorm, mu, X, alpha );
        df = dF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT );
        ratio = f / df;
        X = X - ratio;
        
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

    i_unitVectorMovingFrame = ( rECI_A / norm( rECI_A ));
    k_unitVectorMovingFrame = ( hECI_A / hNorm_A );
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

    
function val = dF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT )

    z = alpha * X^2;
    val = ((r0Norm * v0RadialNorm) / sqrt( mu )) * X^2 * C( z ) + (1 - alpha * r0Norm ) * X^3 * S( z ) + r0Norm * X - sqrt( mu ) * deltaT;
    
end

function val = F( r0Norm, v0RadialNorm, mu, X, alpha )

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

