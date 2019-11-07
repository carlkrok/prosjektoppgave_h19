function [ deltaVStart, deltaVEnd, vIntersectOrbit ] = interceptOrbit( posStart, vStart, posEnd, vEnd, deltaTime, orbitType, mu, tolerance, nMax )


   rStart = norm( posStart );
   rEnd = norm( posEnd );

   deltaTheta = acosd( dot( posStart, posEnd ) / (rStart * rEnd) );

   wOrbit = cross( posStart, posEnd );
   if ((wOrbit(3) < 0 && orbitType == "prograde") || (~(wOrbit(3) < 0 ) && orbitType == "retrograde"))
       deltaTheta = 360 - deltaTheta;
   end

   A = sind( deltaTheta ) * sqrt( (rStart * rEnd) / (1 - cosd( deltaTheta ) ) );

   % TODO: Find first estimate expression
   z0 = -100;
   while F_func( z0, rStart, rEnd, A, mu, deltaTime ) < 0
       z0 = z0 + 0.1;
   end
   z1 = newtonStepZ( z0, rStart, rEnd, A, mu, deltaTime );

   newtonCount = 0;
   while (abs( z1 - z0 ) > tolerance) && (newtonCount < nMax)

       z0 = z1;
       z1 = newtonStepZ( z0, rStart, rEnd, A, mu, deltaTime );
       newtonCount = newtonCount + 1;

   end

   z = z1;


   % Calculating Lagrange constants
   f = 1 - (y_func( z, rStart, rEnd, A ) / rStart);
   g = A * sqrt( y_func( z, rStart, rEnd, A ) / mu);
   df = (sqrt( mu ) / (rStart * rEnd)) * sqrt( y_func( z, rStart, rEnd, A ) / C_stumpff( z ) ) * (z * S_stumpff( z ) - 1);
   dg = 1 - (y_func( z, rStart, rEnd, A ) / rEnd);

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