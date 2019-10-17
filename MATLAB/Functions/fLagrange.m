% Based on eq. 5.31 in [Curtis2011]
function f = fLagrange( r0Norm, alpha, X )

   f = 1 - (X^2 / r0Norm) * C_stumpff( alpha * X^2 );

end