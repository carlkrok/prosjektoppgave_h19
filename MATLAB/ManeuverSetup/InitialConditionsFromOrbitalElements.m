
run earthParametersHPOP;

T = 0; % 
i = 98; % 
O = 0;
w = 0;

rSat_a = rEarth + 632
rSat_p = rEarth + 615

aSat = 0.5 * (rSat_a + rSat_p)

tSat_check = 2*(pi/sqrt(muEarth))*aSat^(3/2)

eSat = 1 - (rSat_p / aSat)

hNorm = sqrt(muEarth * aSat * ( 1 - eSat^2)) 

r0PQW = positionVectorPQW( muEarth, hNorm, eSat, T );
v0PQW = velocityVectorPQW( muEarth, hNorm, eSat, T );
QmatPQWtoECI = transformPQWtoECI( i, O, w );
r0ECI = QmatPQWtoECI * r0PQW
v0ECI = QmatPQWtoECI * v0PQW

rNorm = norm(r0ECI)
vNorm = norm(v0ECI)