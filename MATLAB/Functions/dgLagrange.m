% Based on eq. 5.31 in [Curtis2011]
function dg = dgLagrange( rNorm, alpha, X )

   dg = 1 - (X^2 / rNorm) * C_stumpff( alpha * X^2 );
   %dg = 1 - ( alpha / rNorm) * ( 1 - cos( X / sqrt( alpha ) ) );

end