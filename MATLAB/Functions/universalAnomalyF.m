% Based on eq. 3.59 in [Curtis2011]
function val = universalAnomalyF( r0Norm, v0RadialNorm, mu, X, alpha, deltaT )

   z = alpha * X^2;
   val = ((r0Norm * v0RadialNorm) / sqrt( mu )) * X^2 * C_stumpff( z ) + (1 - alpha * r0Norm ) * X^3 * S_stumpff( z ) + r0Norm * X - sqrt( mu ) * deltaT;

end
