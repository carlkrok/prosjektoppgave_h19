% Based on eq. 3.64 in [Curtis2011]
function val = uninversalAnomalydF( r0Norm, v0RadialNorm, mu, X, alpha )

   z = alpha * X^2;
   val = ((r0Norm * v0RadialNorm) / sqrt( mu )) * X * (1 - alpha * X^2 * S_stumpff( z )) + (1 - alpha * r0Norm) * X^2 * C_stumpff( z ) + r0Norm;

end