function acc = gravityJ2ECI( pos, mu, J_2, R_e )

    r = norm(pos);
    
    acc = - ( mu / (r^3) ) * pos + (( 3 * mu * J_2 * R_e^2 ) / ( 2 * r^5 )) * ( ( 5 * (pos(3))^2 / r^2 ) * [ 1 1 1 ]' - [ 1 1 3 ]' ) * pos
    
end

