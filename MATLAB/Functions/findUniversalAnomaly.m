function X = findUniversalAnomaly( mu, alpha, r0Norm, v0RadialNorm, deltaT, tolerance, nMax )

   X = sqrt( mu ) * abs( alpha ) * deltaT;

   df = uninversalAnomalydF( r0Norm, v0RadialNorm, mu, X, alpha );
   f = universalAnomalyF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT );

   ratio = f / df;

   n = 0;
   while( abs(ratio) > tolerance && n < nMax )

       n = n + 1;
       df = uninversalAnomalydF( r0Norm, v0RadialNorm, mu, X, alpha );
       f = universalAnomalyF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT );
       ratio = f / df;
       X = X - ratio;

   end

   if n == nMax
       %fprintf('\n **No. iterations of Kepler''s equation = %g', n)
       %fprintf('\n   f/df                                = %g\n', f/df)
   end

end