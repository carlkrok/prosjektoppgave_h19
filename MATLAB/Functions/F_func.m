function val = F_func( z, r1, r2, A, mu, deltaT )

   val = (y_func( z, r1, r2, A ) / C_stumpff( z ))^(3 / 2) * S_stumpff( z ) + A * sqrt( y_func( z, r1, r2, A ) ) - sqrt( mu ) * deltaT;

end