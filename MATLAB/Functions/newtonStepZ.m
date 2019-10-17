function z1 = newtonStepZ( z0, r1, r2, A, mu, deltaT )

   z1 = z0 - (F_func( z0, r1, r2, A, mu, deltaT ) / dF_func( z0, r1, r2, A ));

end