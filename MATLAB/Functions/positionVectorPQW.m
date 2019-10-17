function r = positionVectorPQW( mu, h, e, Theta )

   r = ( h^2 / mu ) * ( 1 / ( 1 + e * cosd( Theta ))) * [ cosd( Theta ); sind( Theta ); 0 ];

end

