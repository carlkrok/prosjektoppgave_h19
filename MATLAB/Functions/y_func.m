function val = y_func( z, r1, r2, A )

   val = r1 + r2 + A * ((z * S_stumpff( z ) - 1) / sqrt( C_stumpff( z ) ));
end