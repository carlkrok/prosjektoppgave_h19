function T = orbitPeriod( mu, h, e )

   T  = (( 2 * pi ) / mu^2 ) * (h / sqrt( 1 - e^2 ))^3;

end