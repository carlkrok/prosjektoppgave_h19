function [rLVLH_RelX, rLVLH_RelY, rLVLH_RelZ, rLVLH_RelNorm, sampleT, lastECIPos_B, lastECIVel_B ] = relativeTrajectory( r0_A, v0_A, r0_B, v0_B, anomalyTolerance, nMax, timeHorizon, numPeriods, numSamples, muEarth )

   r_A = r0_A;
   v_A = v0_A;

   r_B = r0_B;
   v_B = v0_B;

   currTime = 0;
   endTime = currTime + timeHorizon; % * numPeriods;
   sampleInterval = (timeHorizon) / numSamples; %  * numPeriods

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

       [r_A, v_A] = nextStateTimeStep( muEarth, r0_A, v0_A, currTime, anomalyTolerance, nMax );
       [r_B, v_B] = nextStateTimeStep( muEarth, r0_B, v0_B, currTime, anomalyTolerance, nMax );

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

