%% LOAD PARAMETERS

% All parameters based on km and s

radiusEarth = 6378; 
gravParamEarth = 398600;

r0Debris = [2000 + radiusEarth, 0, 0];
v0Debris = [0, 8, 0];
hDebris = 0;
iDebris = 0;
raDebris = 0;
eDebris = 0.5;
perDebris = 0;
thetaDebris = 0;


%% (PRELIMINARY) TWO BODY SIMULATION

r0 = r0Debris;
v0 = v0Debris;
h = norm(cross(r0, v0));
mu = gravParamEarth;
e = eDebris;

angleDiff = 1;
numPoints = 360 / angleDiff;

theta = 0;
rEarthX = zeros( numPoints );
rEarthY = zeros( numPoints );

for counter = [1 : numPoints]
    rEarthX( counter ) = radiusEarth * cosd( theta );
    rEarthY( counter ) = radiusEarth * sind( theta );
    theta = theta + angleDiff;
end

figure(1);
hold on;
plot(rEarthX, rEarthY);


for e = [0.1 : 0.05 : 0.4]
    [rX, rY] = orbitPoints(h, mu, e, numPoints);
    plot( rX, rY );
end

hold off;


%% ORBIT DETERMINATION
 
% 1. Calculate r 1 and r 2 using Equation 5.24.
% 2. Choose either a prograde or a retrograde trajectory and calculate ? ? using Equation 5.26.
% 3. Calculate A in Equation 5.35.
% 4. By iteration, using Equations 5.40, 5.43 and 5.45, solve Equation 5.39 for z . The sign of z tells us
% whether the orbit is a hyperbola (z < 0), parabola ( z = 0) or ellipse ( z > 0).
% 5. Calculate y using Equation 5.38.
% 6. Calculate the Lagrange f , g and g functions using Equations 5.46.
% 7. Calculate v 1 and v 2 from Equations 5.28 and 5.29.
% 8. Use r 1 and v 1 (or r 2 and v 2 ) in Algorithm 4.2 to obtain the orbital elements.


startPosition = [5000, 10000, 2100];
endPosition = [-14600, 2500, 7000];
deltaT = 3600;

r1 = norm( startPosition );
r2 = norm( endPosition );


deltaTheta = acosd( dot( startPosition, endPosition ) / (r1 * r2) );

vecNorm = cross( startPosition, endPosition );
% Prograde trajectory
if vecNorm(3) < 0
    % % Retrograde trajectory
    % if ~( cross( r1, r2 )(3) < 0 )
    deltaTheta = 360 - deltaTheta;
end


A = sind( deltaTheta ) * sqrt( (r1 * r2) / (1 - cosd( deltaTheta ) ) );

% TODO: Find first estimate expression
z0 = 1; 
z1 = newtonStep( z0, r1, r2, A, mu, deltaT );

tolerance = 0.00001;
while abs( z1 - z0 ) > tolerance

    z0 = z1;
    z1 = newtonStep( z0, r1, r2, A, mu, deltaT );

end

z = z1;

% Calculating Lagrange constants
f = 1 - (y( z, r1, r2, A ) / r1);
g = A * sqrt( y( z, r1, r2, A ) / mu);
df = (sqrt( mu ) / (r1 * r2)) * sqrt( y( z, r1, r2, A ) / C( z ) ) * (z * S( z ) - 1);
dg = 1 - (y( z, r1, r2, A ) / r2);


% Calculating required velocity in startPosition to arrive at endPosition
% in deltaT seconds
v2 = (1 / g) * (dg * endPosition - startPosition)



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
    

function val = y( z, r1, r2, A )

    val = r1 + r2 + A * ((z * S( z ) - 1) / sqrt( C( z ) ));
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

function [rX, rY] = orbitPoints(h, mu, e, numPoints)

    r = zeros( numPoints );
    rX = zeros( numPoints );
    rY = zeros( numPoints );
    theta = 0;
    angleDiff = 360 / numPoints;
    for counter = [1 : numPoints]
        r( counter ) = ( h^2 / mu ) * 1 / ( 1 + e * cosd( theta ) );
        rX( counter ) = r( counter ) * cosd( theta );
        rY( counter ) = r( counter ) * sind( theta );
        theta = theta + angleDiff;
    end

end
