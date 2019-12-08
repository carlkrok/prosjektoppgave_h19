
function [rXECI, rYECI, rZECI, vXECI, vYECI, vZECI, sampleT] = ECITrajectory( r0, v0, anomalyTolerance, nMax, timeHorizon, numPeriods, numSamples, muEarth )

   if numSamples < 2
       return
   end
   
   currTime = 0;
   endTime = currTime + timeHorizon * numPeriods;
   sampleInterval = (timeHorizon * numPeriods) / (numSamples-1);

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

       [r, v] = nextStateTimeStep( muEarth, r0, v0, currTime, anomalyTolerance, nMax );

       rXECI( sampleIter ) = r( 1 );
       rYECI( sampleIter ) = r( 2 );
       rZECI( sampleIter ) = r( 3 );
       vXECI( sampleIter ) = v( 1 );
       vYECI( sampleIter ) = v( 2 );
       vZECI( sampleIter ) = v( 3 );
       sampleT( sampleIter ) = currTime;

   end

end