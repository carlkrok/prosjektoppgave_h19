function [r, v] = nextStateTimeStep( mu, r0, v0, deltaT, anomalyTolerance, nMax )

   v0Radial = dot( r0, v0 ) / norm( r0 );
   v0RadialNorm = norm( v0Radial );

   r0Norm = norm( r0 );
   v0Norm = norm( v0 );

   alpha = 2 / r0Norm - v0Norm^2 / mu;

   X = findUniversalAnomaly( mu, alpha, r0Norm, v0RadialNorm, deltaT, anomalyTolerance, nMax );

   f = fLagrange( r0Norm, alpha, X );
   g = gLagrange( mu, alpha, X, deltaT );

   r = f * r0 + g * v0;
   rNorm = norm( r );

   df = dfLagrange( r0Norm, rNorm, mu, alpha, X );
   dg = dgLagrange( rNorm, alpha, X );

   v = df * r0 + dg * v0;

end