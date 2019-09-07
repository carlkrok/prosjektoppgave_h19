%% LOAD PARAMETERS

double radiusEarth = 6378; 
double gravParamEarth = 398 600;

double r0Debris = 
double v0Debris = 
double hDebris = 0;
double iDebris = 0;
double raDebris = 0;
double eDebris = 0.2;
double perDebris = 0;
double thetaDebris = 0;

%% ORBIT SIMULATION

double angleDiff = 0.1;

r0 = r0Debris;
v0 = v0Debris;
h = norm(cross(r0, v0));
mu = gravParamEarth;
r = ( h^2 / mu ) * 1 / ( 1 + ( h^2 / ( mu * r0 ) - 1 ) * cosd( 






