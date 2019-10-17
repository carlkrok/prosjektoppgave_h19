function v = velocityVectorPQW( mu, h, e, Theta )

   v = ( mu / h ) * [ -sind( Theta ); e + cosd( Theta ); 0 ];

end

