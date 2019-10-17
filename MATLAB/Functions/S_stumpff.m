function val = S_stumpff( z )

   if z > 0
       val = (sqrt( z ) - sin( sqrt( z ) )) / (sqrt( z ))^3;
   elseif z < 0
       val = (sinh( sqrt( -z ) ) - sqrt( -z )) / (sqrt( -z ))^3;
   else
       val = 1 / 6;
   end

end