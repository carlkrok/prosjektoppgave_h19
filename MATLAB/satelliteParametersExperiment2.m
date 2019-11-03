% Parameters based on Example 7.4 in Curtis

run earthParameters;

%% Chaser

rChaser_a = rEarth + 513.86;
rChaser_p = rEarth + 320.06;
aChaser = 0.5 * (rChaser_a + rChaser_p);
eChaser = 1 - (rChaser_p / aChaser);
hNormChaser = sqrt(muEarth * aChaser * ( 1 - eChaser^2));
tChaser = 1.5484*3600
tChaser_check = 2*(pi/sqrt(muEarth))*aChaser^(3/2)

TChaser = 349.65; % 
iChaser = 40.130; %97.4; % Deg
OChaser = 19.819;
wChaser = 70.662;

rChaserPQW = positionVectorPQW( muEarth, hNormChaser, eChaser, TChaser );
vChaserPQW = velocityVectorPQW( muEarth, hNormChaser, eChaser, TChaser );
QmatPQWtoECIChaser = transformPQWtoECI( iChaser, OChaser, wChaser );
rChaserECI_check = QmatPQWtoECIChaser * rChaserPQW
vChaserECI_check = QmatPQWtoECIChaser * vChaserPQW

hNormChaser_check = norm(rChaserECI_check) * norm(vChaserECI_check);

rChaserECI = [ 1612.75, 5310.19, 3750.33 ]'
vChaserECI = [ -7.35179, 0.463828, 2.46906 ]'

areaSolarChaser = 0.2;
areaDragChaser = 0.1;
CdChaser = 2.2;
CrChaser = 1.0;
massChaser = 2.0;

%% Target

rTarget_a = rEarth + 300;
rTarget_p = rEarth + 300;
aTarget = 0.5 * (rTarget_a + rTarget_p);
eTarget = 1 - (rTarget_p / aTarget);
hNormTarget = sqrt(muEarth * aTarget * ( 1 - eTarget^2));
tTarget = 1.5086*3600
tTarget_check = 2*(pi/sqrt(muEarth))*aTarget^(3/2)

TTarget = 60; % 
iTarget = 40; %97.4; % Deg
OTarget = 20;
wTarget = 0;

rTargetPQW = positionVectorPQW( muEarth, hNormTarget, eTarget, TTarget );
vTargetPQW = velocityVectorPQW( muEarth, hNormTarget, eTarget, TTarget );
QmatPQWtoECITarget = transformPQWtoECI( iTarget, OTarget, wTarget );
rTargetECI_check = QmatPQWtoECITarget * rTargetPQW
vTargetECI_check = QmatPQWtoECITarget * vTargetPQW

hNormTarget_check = norm(rTargetECI_check) * norm(vTargetECI_check);

rTargetECI = [ 1622.39, 5305.10, 3717.55 ]'
vTargetECI = [ -7.29936, 0.492329, 2.48304 ]'










