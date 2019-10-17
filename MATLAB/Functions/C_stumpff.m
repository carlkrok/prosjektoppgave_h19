function val = C_stumpff( z )

   if z > 0
       val = (1 - cos( sqrt( z ) )) / z;
   elseif z < 0
       val = (cosh( sqrt( -z ) ) - 1) / -z;
   else
       val = 1 / 2;
   end

end