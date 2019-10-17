% Program calculates thruster momentum and angle required to reach posEnd
% given posStart, vStart, deltaTime and massSatellite. Assumes thruster
% points in opposite direction of vStart, and operates for tThrust seconds.
% All parameters based on kg, deg, km and s
% Theta defined as angle in XY plane
% Phi defined as angle in ZX plane

run earthParameters;

% For plotting of orbit function
r0Debris = [2000 + rEarth, 0, 0];
v0Debris = [0, 8, 0];
eDebris = 0.2;
angleDiff = 1;

% Lambert's problem, thruster required function
massSatellite = 1.3;
tThrust = 0.1;

posStart = [5000, 10000, 2100];
vStart = [-1.9925, 1.9254, 3.2456];
posEnd = [-14600, 2500, 7000];
vEnd = [-1.9925, 1.9254, 3.2456];

deltaTime = 3600;
orbitTypeDebris = "prograde"; % Alt: retrograde


[ deltaVStart, deltaVEnd ] = changeOrbit( posStart, vStart, posEnd, vEnd, deltaTime, orbitTypeDebris, muEarth )


minTime = 1800;
maxTime = 36000;
deltaTime = 900;

numpoints = ( maxTime - minTime ) / deltaTime;

firstManouver = [ numpoints ];
secondManouver = [ numpoints ];

for manouverTime = minTime : deltaTime : maxTime
    
    [ deltaVStart, deltaVEnd ] = changeOrbit( posStart, vStart, posEnd, vEnd, manouverTime, orbitTypeDebris, muEarth );
    
    firstManouver( ( manouverTime - minTime + deltaTime ) / deltaTime ) = norm( deltaVStart );
    secondManouver( ( manouverTime - minTime + deltaTime ) / deltaTime ) = norm( deltaVEnd );
    
end

figure(3)
hold on
plot(firstManouver)
plot(secondManouver)
legend('First Manouver', 'Second Manouver')
hold off

%plotOrbitWithEarth( r0Debris, v0Debris, eDebris, muEarth, rEarth, angleDiff )
%[ fThruster, thetaThruster, phiThruster ] = thrustRequiredToPoint( posStart, vStart, posEnd, deltaTime, massSatellite, tThrust, orbitTypeDebris, muEarth )




function [ chasersDeltaVStart, chasersDeltaVEnd ] = changeFormation( chasersPosStart, chasersVStart, chasersPosEnd, chasersVEnd, deltaTime, orbitTypeFormation, muEarth )

    numChasers = len( chasersPosStart );
    
    chasersDeltaVStart = [ numChasers ];
    chasersDeltaVEnd = [ numChasers ];
    
    for chaserIterator = 0 : 1 : numChasers
        
        [ chasersDeltaVStart( chaserIterator ), chasersDeltaVEnd( chaserIterator ) ] = changeOrbit( chasersPosStart( chaserIterator ), chasersVStart( chaserIterator ), chasersPosEnd( chaserIterator ), chasersVEnd( chaserIterator ), deltaTime, orbitTypeFormation, muEarth );
        
    end

end

function [ deltaVStart, deltaVEnd ] = changeOrbit( posStart, vStart, posEnd, vEnd, deltaTime, orbitTypeDebris, muEarth )


    rStart = norm( posStart );
    rEnd = norm( posEnd );

    deltaTheta = acosd( dot( posStart, posEnd ) / (rStart * rEnd) );

    wOrbit = cross( posStart, posEnd );
    if ((wOrbit(3) < 0 && orbitTypeDebris == "prograde") || (~(wOrbit(3) < 0 ) && orbitTypeDebris == "retrograde"))
        deltaTheta = 360 - deltaTheta;
    end

    A = sind( deltaTheta ) * sqrt( (rStart * rEnd) / (1 - cosd( deltaTheta ) ) );

    % TODO: Find first estimate expression
    z0 = 0; 
    z1 = newtonStepZ( z0, rStart, rEnd, A, muEarth, deltaTime );
    
    tolerance = 0.00001;
    while abs( z1 - z0 ) > tolerance

        z0 = z1;
        z1 = newtonStepZ( z0, rStart, rEnd, A, muEarth, deltaTime );

    end

    z = z1;


    % Calculating Lagrange constants
    f = 1 - (y_func( z, rStart, rEnd, A ) / rStart);
    g = A * sqrt( y_func( z, rStart, rEnd, A ) / muEarth);
    df = (sqrt( muEarth ) / (rStart * rEnd)) * sqrt( y_func( z, rStart, rEnd, A ) / C_stumpff( z ) ) * (z * S_stumpff( z ) - 1);
    dg = 1 - (y_func( z, rStart, rEnd, A ) / rEnd);

    % Calculating required velocity in startPosition to arrive at end orbit
    % in deltaT seconds
    vRequiredStart = (1 / g) * (posEnd - f * posStart);
    vIntersectOrbit = (1 / g) * ( dg * posEnd - posStart );

    % Change in velocity required
    deltaVStart = vStart - vRequiredStart;
    deltaVEnd = vIntersectOrbit - vEnd;


end

function [ fThruster, thetaThruster, phiThruster ] = thrustRequiredToPoint( posStart, vStart, posEnd, deltaTime, massSatellite, tThrust, orbitTypeDebris, muEarth )

    rStart = norm( posStart );
    rEnd = norm( posEnd );

    deltaTheta = acosd( dot( posStart, posEnd ) / (rStart * rEnd) );

    wOrbit = cross( posStart, posEnd );
    if ((wOrbit(3) < 0 && orbitTypeDebris == "prograde") || (~(wOrbit(3) < 0 ) && orbitTypeDebris == "retrograde"))
        deltaTheta = 360 - deltaTheta;
    end

    A = sind( deltaTheta ) * sqrt( (rStart * rEnd) / (1 - cosd( deltaTheta ) ) );

    % TODO: Find first estimate expression
    z0 = 0; 
    z1 = newtonStepZ( z0, rStart, rEnd, A, muEarth, deltaTime );

    tolerance = 0.00001;
    while abs( z1 - z0 ) > tolerance

        z0 = z1;
        z1 = newtonStepZ( z0, rStart, rEnd, A, muEarth, deltaTime );

    end

    z = z1;


    % Calculating Lagrange constants
    f = 1 - (y_func( z, rStart, rEnd, A ) / rStart);
    g = A * sqrt( y_func( z, rStart, rEnd, A ) / muEarth);
    df = (sqrt( muEarth ) / (rStart * rEnd)) * sqrt( y_func( z, rStart, rEnd, A ) / C_stumpff( z ) ) * (z * S_stumpff( z ) - 1);
    dg = 1 - (y_func( z, rStart, rEnd, A ) / rEnd);

    % Calculating required velocity in startPosition to arrive at endPosition
    % in deltaT seconds
    vRequired = (1 / g) * (posEnd - f * posStart);

    % Change in velocity required
    deltaV = vStart - vRequired;

    % Momentum that thruster must provide
    deltaH = deltaV * massSatellite;

    fThruster = norm( deltaH ) / tThrust; % Newtons

    % Theta defined as angle in XY plane
    thetaThruster = acosd( dot( vStart(1:2) , vRequired(1:2) ) / (norm( vStart(1:2) ) * norm( vRequired(1:2) )));

    % Phi defined as angle in ZX plane
    phiThruster = acosd( dot( vStart(1,3) , vRequired(1,3) ) / (norm( vStart(1,3) ) * norm( vRequired(1,3) )));

end

function plotOrbitWithEarth( r0, v0, e, muEarth, rEarth, angleDiff )

    h = norm(cross(r0, v0));
    numPoints = 360 / angleDiff;
    theta = 0;
    rEarthX = zeros( numPoints );
    rEarthY = zeros( numPoints );

    for counter = [1 : numPoints]
        rEarthX( counter ) = rEarth * cosd( theta );
        rEarthY( counter ) = rEarth * sind( theta );
        theta = theta + angleDiff;
    end

    figure(1);
    hold on;
    plot(r0( 1 ), r0( 2 ), '*');
    plot(rEarthX, rEarthY);

    [rX, rY] = orbitPoints(h, muEarth, e, numPoints);
    
    plot( rX, rY );
    hold off;
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


% TODOs
% add 'selection of the least costly trajectory that satisfies also the constraints'. Cost indexes = minimum time, minimum Delta v, do a convex combination
% add the computation of the propagation of the uncertainties (from the actuation being different from the thought one to environmental phenomena, like J2, atmospheric drag, etc to distribution of the final orbital parameters, thrust direction, thrust performance, uncertainty on deltaV)
% CW equations
% Error parameters - max & min - simulate and see corresponding output
% error, assuming Keplerian motion




