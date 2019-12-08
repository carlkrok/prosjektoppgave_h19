% Based on eq. 5.31 in [Curtis2011]
function df = dfLagrange( r0Norm, rNorm, mu, alpha, X )

   df = sqrt( mu ) / (r0Norm * rNorm) * (alpha * X^3 * S_stumpff( alpha * X^2 ) - X) ;
   %df =  - (sqrt( mu * alpha ) / (r0Norm * rNorm)) * (alpha * X^3 * sin( X / sqrt( alpha ) )) ;

end