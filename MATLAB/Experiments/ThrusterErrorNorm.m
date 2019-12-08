
error10 = 1.89012105481708;
error8 = 1.20968626290305;
error6 = 0.680454032785321;
error4 = 0.302427557186588;
error2 = 0.0756094925139242;

figure(1)
hold on
grid on
title('Norm of End Position Error')
plot( [10, 8, 6, 4, 2],[error10, error8, error6, error4, error2] )
legend('HPOP Keplerian Orbit')
xlabel('Thrust Time [s]')
ylabel('Distance [m]')
hold off