% Based on eq. 5.31 in [Curtis2011]
function g = gLagrange( mu, alpha, X, deltaT )

   g = deltaT - (1 / sqrt( mu )) * X^3 * S_stumpff( alpha * X^2 );
   %g = deltaT - ( alpha / sqrt( mu )) * ( X *  - sqrt( alpha ) * sin( X / sqrt( alpha ) ) );

end
