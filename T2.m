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
i_B = 50;
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


hECI_A = cross( rECI_A, vECI_A );

i_unitVectorMovingFrame = ( rECI_A / norm( rECI_A ));
k_unitVectorMovingFrame = ( hECI_A / hNorm_A );
j_unitVectorMovingFrame = cross( k_unitVectorMovingFrame, i_unitVectorMovingFrame );


QmatECItoLVLH = [i_unitVectorMovingFrame';
                 j_unitVectorMovingFrame';
                 k_unitVectorMovingFrame']
                    
angularVelocityLVLH = hECI_A / (norm( rECI_A ))^2
angularAccelerationLVLH = -2 * (dot( vECI_A, rECI_A ) / (norm( rECI_A ))^2) * angularVelocityLVLH

aECI_A = -(muEarth / norm(rECI_A)^3) * rECI_A
aECI_B = -(muEarth / norm(rECI_B)^3) * rECI_B


rECI_Rel = rECI_B - rECI_A;
vECI_Rel = vECI_B - vECI_A - cross( angularVelocityLVLH, rECI_Rel );
aECI_Rel = aECI_B - aECI_A - cross( angularAccelerationLVLH, rECI_Rel ) - cross( angularVelocityLVLH, cross( angularVelocityLVLH, rECI_Rel )) - 2 * cross( angularVelocityLVLH, vECI_Rel );


rLVLH_Rel = QmatECItoLVLH * rECI_Rel;
vLVLH_Rel = QmatECItoLVLH * vECI_Rel;
aLVLH_Rel = QmatECItoLVLH * aECI_Rel;


%% Plot of relative motion

orbitPeriod_A = orbitPeriod( muEarth, hNorm_A, e_A )


numSamples = 1000;
sampleInterval = orbitPeriod / numSamples;




function X = findAnomaly( mu, alpha, r0, v0_radial, deltaT )

    X0 = sqrt( mu ) * abs( alpha ) * deltaT;
    
    

end


    
function val = dF( )

    
end

function val = F( )


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


function [r, v] = nextState( mu, r0, v0, deltaT )

    v0_radial = dot( r0, v0 ) / norm( r0 );
    
    alpha = 2 / norm( r0 ) - norm( v0 )^2 / mu;
    
    

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







