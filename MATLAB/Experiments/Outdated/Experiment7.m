%% Clear Workspace
clc
clear all
format long g
close all


%% Load Parameters

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC 

run earthParametersHPOP; 

r0ECI_target = [ 6990.1366, 0.0, 0.0 ]';
v0ECI_target = [ 0.0, -1.04884571201957, 7.50091780930118 ]';

r0ECI_chaser = [ 6952.1366, 0.0, 0.0 ]';
v0ECI_chaser = [ 0.0, -0.976461903728063, 7.51834019814233 ]';

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0,...
                  'Thrust',0, 'stepCounter', 0, 'velocityChangeLVLH', 0, 'thrustLVLHAcceleration', 0);

anomalyErrorTolerance = 10^(-12);
anomalyMaxIterations = 2000;

SAT_Const
constants
load DE436Coeff.mat
PC = DE436Coeff; % JPL Planetary Ephemeris DECEMBER 14 1949 - DECEMBER 21 2149


% read Earth gravity field coefficients
Cnm = zeros(181,181); % Geopotential coefficient (cos)
Snm = zeros(181,181); % Geopotential coefficient (sin)
fid = fopen('GGM03S.txt','r'); % GRACE Gravity Model 03 (GGM03), GGM03S - complete to harmonic degree 180
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

% read Earth orientation parameters 1999 - 2018
fid = fopen('eop19990101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

% read space weather data 1999 - 2017
fid = fopen('sw19990101.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
swdata = fscanf(fid,'%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4f %2i %4i %6f %2i %6f %6f %6f %6f %6f',[33 inf]);
fclose(fid);

% read space weather data 1997 - 2016 Solar data
fid = fopen('SOLFSMY.txt','r'); % SOLFSMY_2019
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% READ GEOMAGNETIC STORM DTC VALUE
fid = fopen('DTCFILE.txt','r'); % DTCFILE_2019
%  ------------------------------------------------------------------------
% | DTC YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);

% read space weather data
fid = fopen('SOLRESAP.txt','r'); % SOLRESAP_2019
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[12 inf]);
fclose(fid);


% epoch state
year = 2002;
mon = 04;
day = 24;
hour = 00;
min = 00;
sec = 00;
Y0 = [ r0ECI_chaser', v0ECI_chaser'].*10^3;
AuxParam.area_solar = 0.2;
AuxParam.area_drag = 0.1;
AuxParam.mass = 2.0;
AuxParam.Cr = 1.0;
AuxParam.Cd = 2.2;

% epoch
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
%Y0 = ECEF2ECI(Mjd_UTC, Y0); Already in ECI

AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 0;
AuxParam.m       = 0;
AuxParam.sun     = 0;
AuxParam.moon    = 0;
AuxParam.planets = 0;
AuxParam.sRad    = 0;
AuxParam.drag    = 0;
AuxParam.SolidEarthTides = 0;
AuxParam.OceanTides = 0;
AuxParam.Relativity = 0;


%% Experiment Setup

maneuverEndTime = 2500;
maneuverStartTime = 300;
Step   = 0.05;   % [s]
N_Step = round(maneuverEndTime*1/Step); 
N_Step_Initial = round(maneuverStartTime *1/Step);
              
Mjd0   = Mjd_UTC;

% shorten PC, eopdata, swdata, Cnm, and Snm
num = fix(N_Step*Step/86400)+2;
JD = Mjd_UTC+2400000.5;
i = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
PC = PC(i:i+num,:);
mjd = (floor(Mjd_UTC));
i = find(mjd==eopdata(4,:),1,'first');
eopdata = eopdata(:,i:i+num);
i = find((year==swdata(1,:)) & (mon==swdata(2,:)) & (day==swdata(3,:)),1,'first');
swdata = swdata(:,i-3:i+num);
%Cnm = Cnm(1:AuxParam.n+1,1:AuxParam.n+1);
%Snm = Snm(1:AuxParam.n+1,1:AuxParam.n+1);

orbitType_chaser = "retrograde";

numSamplePointsInitialTrajectorySimple = 200;
numSamplePointsFinalTrajectorySimple = 1000;

%% Monte Carlo Experiment Setup

MCsampleNum = 11;

% meanDeviationTimeSetup = 0;
maxDeviationTimeSetup = 3;
% MCdeviationTimes = meanDeviationTimeSetup - maxDeviationTimeSetup + ( 2 * maxDeviationTimeSetup * rand( MCsampleNum, 1 ) );

MCdeviationTimes = (-maxDeviationTimeSetup : (2*maxDeviationTimeSetup)/(MCsampleNum-1) : maxDeviationTimeSetup)';

experimentStartTimeIdeal = maneuverStartTime;


%% Initial Orbit Determination

orbitPeriod_chaser = 3000;

targetY0 = [ r0ECI_target', v0ECI_target' ].*10^3;
[tagetEph] = ephemeris_v3(targetY0, N_Step, Step);
rECIManouverEnd_target = [ tagetEph( N_Step+1, 2 ); tagetEph( N_Step+1, 3 ); tagetEph( N_Step+1, 4 ) ]./10^3;
vECIManouverEnd_target = [ tagetEph( N_Step+1, 5 ); tagetEph( N_Step+1, 6 ); tagetEph( N_Step+1, 7 ) ]./10^3;

% B's ideal position at experiment start
chaserY0 = [ r0ECI_chaser', v0ECI_chaser' ].*10^3;
[chaserEph] = ephemeris_v3(chaserY0, N_Step_Initial, Step);
rECIManeuverStartIdeal_chaser = [ chaserEph( N_Step_Initial+1, 2 ); chaserEph( N_Step_Initial+1, 3 ); chaserEph( N_Step_Initial+1, 4 ) ]./10^3;
vECIManeuverStartIdeal_chaser = [ chaserEph( N_Step_Initial+1, 5 ); chaserEph( N_Step_Initial+1, 6 ); chaserEph( N_Step_Initial+1, 7 ) ]./10^3;

% Required velocity change satellite B
[ deltaVStartECI_chaser, deltaVEndECI_chaser, vIntersectOrbit_chaser ] = interceptOrbit( rECIManeuverStartIdeal_chaser, vECIManeuverStartIdeal_chaser, rECIManouverEnd_target, vECIManouverEnd_target, maneuverEndTime - maneuverStartTime, orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );

[ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartIdeal_chaser, vECIManeuverStartIdeal_chaser );
deltaVStartLVLH_chaser = QmatECItoLVLH_chaser * deltaVStartECI_chaser;

AuxParam.velocityChangeLVLH = deltaVStartLVLH_chaser*1000;


%% Monte Carlo Experimetns


%%%%%%%%%%%%%  SIMPLE MODEL


MCSimplePosEnd = zeros( MCsampleNum, 3);
MCSimplePosStart = zeros( MCsampleNum, 3);

MCSimpleVelEnd = zeros(MCsampleNum, 3);
MCSimpleVelStart = zeros(MCsampleNum, 3);

MCSimpleECI_X_Trajectories = zeros(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);
MCSimpleECI_Y_Trajectories = zeros(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);
MCSimpleECI_Z_Trajectories = zeros(MCsampleNum, numSamplePointsInitialTrajectorySimple + numSamplePointsFinalTrajectorySimple);

% B's position at first iteration start
[ rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + MCdeviationTimes( 1 ), anomalyErrorTolerance, anomalyMaxIterations );

[rXECITrajectoryExperimentStart_chaser, rYECITrajectoryExperimentStart_chaser, rZECITrajectoryExperimentStart_chaser, vXECITrajectoryExperimentStart_chaser, vYECITrajectoryExperimentStart_chaser, vZECITrajectoryExperimentStart_chaser, sampleTECITrajectoryExperimentStart_chaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartTime + MCdeviationTimes( 1 ), 1, numSamplePointsInitialTrajectorySimple, muEarth );
MCSimpleECI_X_Trajectories(1, 1:numSamplePointsInitialTrajectorySimple) = rXECITrajectoryExperimentStart_chaser;
MCSimpleECI_Y_Trajectories(1, 1:numSamplePointsInitialTrajectorySimple) = rYECITrajectoryExperimentStart_chaser;
MCSimpleECI_Z_Trajectories(1, 1:numSamplePointsInitialTrajectorySimple) = rZECITrajectoryExperimentStart_chaser;

MCSimplePosStart( 1, : ) = rECIExperimentStartSimple_chaser';
MCSimpleVelStart( 1, : ) = vECIExperimentStartSimple_chaser';

[ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser );
QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
deltaVExperimentStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;

% B's position at first iteration end
[ rECIExperimentEndSimple_chaser, vECIExperimentEndSimple_chaser ] = nextStateTimeStep( muEarth, rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser + deltaVExperimentStart_chaser, maneuverEndTime - maneuverStartTime - MCdeviationTimes( 1 ), anomalyErrorTolerance, anomalyMaxIterations );

[rXECITrajectoryExperimentEnd_chaser, rYECITrajectoryExperimentEnd_chaser, rZECITrajectoryExperimentEnd_chaser, vXECITrajectoryExperimentEnd_chaser, vYECITrajectoryExperimentEnd_chaser, vZECITrajectoryExperimentEnd_chaser, sampleTECITrajectoryExperimentEnd_chaser] = ECITrajectory( rECIExperimentStartSimple_chaser, vECIExperimentStartSimple_chaser + deltaVExperimentStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverEndTime - maneuverStartTime - MCdeviationTimes( 1 ), 1, numSamplePointsFinalTrajectorySimple, muEarth );
MCSimpleECI_X_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rXECITrajectoryExperimentEnd_chaser;
MCSimpleECI_Y_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rYECITrajectoryExperimentEnd_chaser;
MCSimpleECI_Z_Trajectories(1, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rZECITrajectoryExperimentEnd_chaser;

MCSimplePosEnd( 1, : ) = rECIExperimentEndSimple_chaser';
MCSimpleVelEnd( 1, : ) = vECIExperimentEndSimple_chaser';



for experimentIndex = 2 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );
  
    %%%%%%%%%%%%%  SIMPLE MODEL
    % B's position at maneuver start
    [ rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser ] = nextStateTimeStep( muEarth, r0ECI_chaser, v0ECI_chaser, experimentStartTimeIdeal + thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    MCSimplePosStart( experimentIndex, : ) = rECIManeuverStartSimple_chaser';
    MCSimpleVelStart( experimentIndex, : ) = vECIManeuverStartSimple_chaser';
    
    [rXECITrajectoryManeuverStart_chaser, rYECITrajectoryManeuverStart_chaser, rZECITrajectoryManeuverStart_chaser, vXECITrajectoryManeuverStart_chaser, vYECITrajectoryManeuverStart_chaser, vZECITrajectoryManeuverStart_chaser, sampleTECITrajectoryManeuverStart_chaser] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverStartTime + MCdeviationTimes( experimentIndex ), 1, numSamplePointsInitialTrajectorySimple, muEarth );
    MCSimpleECI_X_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rXECITrajectoryManeuverStart_chaser;
    MCSimpleECI_Y_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rYECITrajectoryManeuverStart_chaser;
    MCSimpleECI_Z_Trajectories(experimentIndex, 1:numSamplePointsInitialTrajectorySimple) = rZECITrajectoryManeuverStart_chaser;

    [ QmatECItoLVLH_chaser ] = ECIToLVLH( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser );
    QmatLVLHtoECI_chaser = QmatECItoLVLH_chaser';
    deltaVManeuverStart_chaser = QmatLVLHtoECI_chaser * deltaVStartLVLH_chaser;

    % B's position at experiment end
    [ rECIExperimentEndSimple_chaser, vECIExperimentEndSimple_chaser ] = nextStateTimeStep( muEarth, rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser + deltaVManeuverStart_chaser, maneuverEndTime - maneuverStartTime - thisDeviationTime, anomalyErrorTolerance, anomalyMaxIterations );

    MCSimplePosEnd( experimentIndex, : ) = rECIExperimentEndSimple_chaser';
    MCSimpleVelEnd( experimentIndex, : ) = vECIExperimentEndSimple_chaser';

    [rXECITrajectoryManeuverEnd_chaser, rYECITrajectoryManeuverEnd_chaser, rZECITrajectoryManeuverEnd_chaser, vXECITrajectoryManeuverEnd_chaser, vYECITrajectoryManeuverEnd_chaser, vZECITrajectoryManeuverEnd_chaser, sampleTECITrajectoryManeuverEnd_chaser] = ECITrajectory( rECIManeuverStartSimple_chaser, vECIManeuverStartSimple_chaser + deltaVManeuverStart_chaser, anomalyErrorTolerance, anomalyMaxIterations, maneuverEndTime - maneuverStartTime - MCdeviationTimes( experimentIndex ), 1, numSamplePointsFinalTrajectorySimple, muEarth );
    MCSimpleECI_X_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rXECITrajectoryManeuverEnd_chaser;
    MCSimpleECI_Y_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rYECITrajectoryManeuverEnd_chaser;
    MCSimpleECI_Z_Trajectories(experimentIndex, (numSamplePointsInitialTrajectorySimple+1):(numSamplePointsInitialTrajectorySimple+numSamplePointsFinalTrajectorySimple)) = rZECITrajectoryManeuverEnd_chaser;
    
end


%%%%%%%%%%%%%  HPOP MODEL


%%%%%%%  Experiment 1 Only Thrust effect d = 1s


AuxParam.thrustDuration = 1;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 0;
AuxParam.m       = 0;
AuxParam.sun     = 0;
AuxParam.moon    = 0;
AuxParam.planets = 0;
AuxParam.sRad    = 0;
AuxParam.drag    = 0;
AuxParam.SolidEarthTides = 0;
AuxParam.OceanTides = 0;
AuxParam.Relativity = 0;


MC_1_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_1_HPOP_PosStart = zeros( MCsampleNum, 3);

MC_1_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_1_HPOP_VelStart = zeros(MCsampleNum, 3);

MC_1_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+3);
MC_1_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+3);
MC_1_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+3);


for experimentIndex = 1 : MCsampleNum
    
    thisDeviationTime = MCdeviationTimes( experimentIndex );

    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    
    N_Step_Initial = floor((maneuverStartTime + thisDeviationTime) *1/Step); % 
    N_Step_Thrust = round(AuxParam.thrustDuration*1/Step); %  
    N_Step_Final = ceil(maneuverEndTime*1/Step) - N_Step_Initial - N_Step_Thrust; %  

    AuxParam.thrustLVLHAcceleration = AuxParam.velocityChangeLVLH./AuxParam.thrustDuration;

    AuxParam.Thrust = 0;
    [initialEph] = ephemeris_v3(Y0, N_Step_Initial, Step);
    currY = [ initialEph(N_Step_Initial+1, 2), initialEph(N_Step_Initial+1, 3), initialEph(N_Step_Initial+1, 4), initialEph(N_Step_Initial+1, 5), initialEph(N_Step_Initial+1, 6), initialEph(N_Step_Initial+1, 7) ];

    AuxParam.Thrust = 1;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime  + MCdeviationTimes( 1 ))/const.DAYSEC);
    [thrustEph] = ephemeris_v3(currY, N_Step_Thrust, Step);
    currY = [ thrustEph(N_Step_Thrust+1, 2), thrustEph(N_Step_Thrust+1, 3), thrustEph(N_Step_Thrust+1, 4), thrustEph(N_Step_Thrust+1, 5), thrustEph(N_Step_Thrust+1, 6), thrustEph(N_Step_Thrust+1, 7) ];

    AuxParam.Thrust = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + AuxParam.thrustDuration + MCdeviationTimes( 1 ))/const.DAYSEC);
    [finalEph] = ephemeris_v3(currY, N_Step_Final, Step);

    Eph = [initialEph; thrustEph; finalEph];

    MC_1_HPOP_PosEnd( experimentIndex, : ) = [Eph(N_Step+3, 2), Eph(N_Step+3, 3), Eph(N_Step+3, 4)]./10^3;
    MC_1_HPOP_VelEnd( experimentIndex, : ) = [Eph(N_Step+3, 5), Eph(N_Step+3, 6), Eph(N_Step+3, 7)]./10^3;

    MC_1_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph(:, 2)'./10^3;
    MC_1_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph(:, 3)'./10^3;
    MC_1_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph(:, 4)'./10^3;
    
end


%% Statistical Analysis

[ QmatECItoLVLH_target ] = ECIToLVLH( rECIManouverEnd_target, vECIManouverEnd_target );

absDeviationEndPosSimple = zeros( MCsampleNum, 1 );
relEndPosECISimple_chaser = MCSimplePosEnd - ones(MCsampleNum,1)*rECIManouverEnd_target';
relEndPosLVLHSimple_chaser = zeros(MCsampleNum,3)
for posIter = 1:MCsampleNum
    thisPos = relEndPosECISimple_chaser( posIter, : )';
    thisPos = QmatECItoLVLH_target * thisPos;
    relEndPosLVLHSimple_chaser( posIter, : ) = thisPos';
end


absDeviationEndPosHPOP_1 = zeros( MCsampleNum, 1 );
relEndPosECIHPOP_1_chaser = MC_1_HPOP_PosEnd - ones(MCsampleNum,1)*rECIManouverEnd_target';
relEndPosLVLHHPOP_1_chaser = zeros(MCsampleNum,3)
for posIter = 1:MCsampleNum
    thisPos = relEndPosECIHPOP_1_chaser( posIter, : )';
    thisPos = QmatECItoLVLH_target * thisPos;
    relEndPosLVLHHPOP_1_chaser( posIter, : ) = thisPos';
end

for dataIndex = 1 : MCsampleNum
    
    thisEndPosSimple = MCSimplePosEnd( dataIndex, : );
    
    thisEndPosHPOP_1 = MC_1_HPOP_PosEnd( dataIndex, : );

    absDeviationEndPosSimple( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosSimple' );
    absDeviationEndPosHPOP_1( dataIndex ) = norm( rECIManouverEnd_target - thisEndPosHPOP_1' );

end


%% Additoinal Trajectories


% A's trajectory
[rXECITrajectory_target, rYECITrajectory_target, rZECITrajectory_target, vXECITrajectory_target, vYECITrajectory_target, vZECITrajectory_target, sampleTECITrajectory_target] = ECITrajectory( r0ECI_target, v0ECI_target, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_chaser, 1, 1000, muEarth );

% B's old trajectory
[rXECITrajectory_chaserOld, rYECITrajectory_chaserOld, rZECITrajectory_chaserOld, vXECITrajectory_chaserOld, vYECITrajectory_chaserOld, vZECITrajectory_chaserOld, sampleTECITrajectory_chaserOld] = ECITrajectory( r0ECI_chaser, v0ECI_chaser, anomalyErrorTolerance, anomalyMaxIterations, orbitPeriod_chaser, 1, 1000, muEarth );



%% Plots


figure(1)
hold on
grid on
title('Error in ECI XYZ-Position of MC Simulations')
plot3( relEndPosECISimple_chaser( :, 1 ), relEndPosECISimple_chaser( :, 2 ), relEndPosECISimple_chaser( :, 3 ), '*' )
plot3( relEndPosECIHPOP_1_chaser( :, 1 ), relEndPosECIHPOP_1_chaser( :, 2 ), relEndPosECIHPOP_1_chaser( :, 3 ), '*' )
plot3(0,0,0,'m+', 'linewidth',8)
legend('Simple Model', 'HPOP 1',  'Goal Position')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
hold off


figure(2)
hold on
grid on
title('Maneuver End Position LVLH')
plot3( relEndPosLVLHSimple_chaser( :, 1 ), relEndPosLVLHSimple_chaser( :, 2 ), relEndPosLVLHSimple_chaser( :, 3 ), '*' )
plot3( relEndPosLVLHHPOP_1_chaser( :, 1 ), relEndPosLVLHHPOP_1_chaser( :, 2 ), relEndPosLVLHHPOP_1_chaser( :, 3 ), '*' )
plot3(0,0,0,'m+', 'linewidth',8)
legend('Simple Model', 'HPOP 1',  'Goal Position')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
hold off


figure(3)
hold on
grid on
title('Norm of Error in Point of Rendezvous of MC Simulations')
plot( absDeviationEndPosSimple)
plot( absDeviationEndPosHPOP_1 )
legend('Simple Model', 'HPOP 1')
xlabel('Sample')
ylabel('Distance [km]')
hold off


figure(4)
hold on
grid on
title('ECI Trajectories')
%axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
%h1 = plot3( rXECITrajectory_target, rYECITrajectory_target, rZECITrajectory_target, '-.')
%h2 = plot3( rXECITrajectory_chaserOld, rYECITrajectory_chaserOld, rZECITrajectory_chaserOld, '-.')
%[ sx, sy, sz ] = sphere;
%surf( sx*rEarth, sy*rEarth, sz*rEarth, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.05 );
for plotIndex = 1 : MCsampleNum
    plot3( MCSimpleECI_X_Trajectories(plotIndex, :), MCSimpleECI_Y_Trajectories(plotIndex, :), MCSimpleECI_Z_Trajectories(plotIndex, :), 'b')
    plot3( MC_1_HPOP_ECI_X_Trajectories(plotIndex, :), MC_1_HPOP_ECI_Y_Trajectories(plotIndex, :), MC_1_HPOP_ECI_Z_Trajectories(plotIndex, :), 'r')
end
%plot3( r0ECI_chaser(1), r0ECI_chaser(2), r0ECI_chaser(3), '*k' )
%plot3( rECIExperimentStartSimple_chaser(1), rECIExperimentStartSimple_chaser(2), rECIExperimentStartSimple_chaser(3), '*k' )
plot3( rECIManouverEnd_target(1), rECIManouverEnd_target(2), rECIManouverEnd_target(3), '*k' )
%text( r0ECI_chaser(1), r0ECI_chaser(2), r0ECI_chaser(3), 'Chaser Experiment Start' )
%text( rECIExperimentStartSimple_chaser(1), rECIExperimentStartSimple_chaser(2), rECIExperimentStartSimple_chaser(3), 'Chaser Ideal Maneuver Start' )
text( rECIManouverEnd_target(1), rECIManouverEnd_target(2), rECIManouverEnd_target(3), 'Chaser Ideal Maneuver End' )
%legend('Simple Model initial point', 'Simple Model maneuver start', 'Simple Model maneuver end')
%legend([h1, h2], 'Target Trajectory', 'Chaser Initial Trajectory')
hold off

