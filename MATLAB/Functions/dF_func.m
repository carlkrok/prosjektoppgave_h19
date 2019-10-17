function val = dF_func( z, r1, r2, A )

   if z == 0
       val = (sqrt( 2 ) / 40) * y_func( 0, r1, r2, A )^(3 / 2) + (A / 8) * (sqrt( y_func( 0, r1, r2, A ) ) + A * sqrt( 1 / (2 * y_func( 0, r1, r2, A ) ) ));
   else
       val = (y_func( z, r1, r2, A ) / C_stumpff( z ))^(3/2) * ((1 / (2 * z)) * (C_stumpff( z ) - (3 / 2) * (S_stumpff( z )^2 / C_stumpff( z ))) + (3 / 4) * (S_stumpff( z )^2 / C_stumpff( z ))) + (A / 8) * (3 * (S_stumpff(z ) / C_stumpff( z )) * sqrt( y_func( z, r1, r2, A ) ) + A * sqrt( C_stumpff( z ) / y_func( z, r1, r2, A )));
   end

end