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


%% Experiment Setup

maneuverEndTime = 2500;
maneuverStartTime = 500;
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

MCsampleNum = 3000;

meanDeviationTimeSetup = 0;
%maxDeviationTimeSetup = 1.5;
stdDeviationTimeSetup = 0.2;

%MCtimeDeviation = meanDeviationTimeSetup - stdDeviationTimeSetup + ( 2 * stdDeviationTimeSetup * rand( MCsampleNum, 1 ) );
%MCtimeDeviation = (meanDeviationTimeSetup-maxDeviationTimeSetup : (2*maxDeviationTimeSetup)/(MCsampleNum-1) : meanDeviationTimeSetup+maxDeviationTimeSetup)';
MCtimeDeviation = meanDeviationTimeSetup + stdDeviationTimeSetup .* randn( MCsampleNum, 1 );


meanThrustOutputFactor = 1;
stdThrustOutputUncertaintyFactor = 0.0025;

%MCthrustOutputDeviation = meanThrustOutputFactor - stdThrustOutputUncertaintyFactor + ( 2 * stdThrustOutputUncertaintyFactor * rand( MCsampleNum, 1 ) );
%MCthrustOutputDeviation = (meanThrustOutputFactor-stdThrustOutputUncertaintyFactor : (2*stdThrustOutputUncertaintyFactor)/(MCsampleNum-1) : meanThrustOutputFactor+stdThrustOutputUncertaintyFactor)';
MCthrustOutputDeviation = meanThrustOutputFactor + stdThrustOutputUncertaintyFactor .* randn( MCsampleNum, 1 );


meanThrustDirectionErrorRollDeg = 0;
stdThrustDirectionErrorRollDeg = 0.001;

%MCthrustDirectionRollDeviationDeg = meanThrustDirectionErrorRollDeg - stdThrustDirectionErrorRollDeg + ( 2 * stdThrustDirectionErrorRollDeg * rand( MCsampleNum, 1 ) );
%MCthrustDirectionRollDeviationRad = (meanThrustDirectionErrorRollDeg-maxThrustDirectionErrorRollDeg : (2*maxThrustDirectionErrorRollDeg)/(MCsampleNum-1) : meanThrustDirectionErrorRollDeg+maxThrustDirectionErrorRollDeg)'*pi/180;
MCthrustDirectionRollDeviationDeg = meanThrustDirectionErrorRollDeg + stdThrustDirectionErrorRollDeg .* randn( MCsampleNum, 1 );


meanThrustDirectionErrorPitchDeg = 0;
stdThrustDirectionErrorPitchDeg = 0.001;

%MCthrustDirectionPitchDeviationDeg = meanThrustDirectionErrorPitchDeg - stdThrustDirectionErrorPitchDeg + ( 2 * stdThrustDirectionErrorPitchDeg * rand( MCsampleNum, 1 ) );
%MCthrustDirectionPitchDeviationRad = (meanThrustDirectionErrorPitchDeg-maxThrustDirectionErrorPitchDeg : (2*maxThrustDirectionErrorPitchDeg)/(MCsampleNum-1) : meanThrustDirectionErrorPitchDeg+maxThrustDirectionErrorPitchDeg)'*pi/180;
MCthrustDirectionPitchDeviationDeg = meanThrustDirectionErrorPitchDeg + stdThrustDirectionErrorPitchDeg .* randn( MCsampleNum, 1 );


meanThrustDirectionErrorYawDeg = 0;
stdThrustDirectionErrorYawDeg = 0.001;

%MCthrustDirectionYawDeviationDeg = meanThrustDirectionErrorYawDeg - stdThrustDirectionErrorYawDeg + ( 2 * stdThrustDirectionErrorYawDeg * rand( MCsampleNum, 1 ) );
%MCthrustDirectionYawDeviationRad = (meanThrustDirectionErrorYawDeg-maxThrustDirectionErrorYawDeg : (2*maxThrustDirectionErrorYawDeg)/(MCsampleNum-1) : meanThrustDirectionErrorYawDeg+maxThrustDirectionErrorYawDeg)'*pi/180;
MCthrustDirectionYawDeviationDeg = meanThrustDirectionErrorYawDeg + stdThrustDirectionErrorYawDeg .* randn( MCsampleNum, 1 );



%% Initial Orbit Determination

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

targetY0 = [ r0ECI_target', v0ECI_target' ].*10^3;
targetEph = ephemeris_v3(targetY0, N_Step, Step);
rECIManouverEnd_target = [ targetEph( N_Step+1, 2 ); targetEph( N_Step+1, 3 ); targetEph( N_Step+1, 4 ) ]./10^3;
vECIManouverEnd_target = [ targetEph( N_Step+1, 5 ); targetEph( N_Step+1, 6 ); targetEph( N_Step+1, 7 ) ]./10^3;

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


%% Experiments

AuxParam.n       = 40;
AuxParam.m       = 40;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.sRad    = 1;
AuxParam.drag    = 1;
AuxParam.SolidEarthTides = 0;
AuxParam.OceanTides = 0;
AuxParam.Relativity = 0;

targetPreciseEph = ephemeris_v3(targetY0, N_Step, Step);
rECIManouverEnd_targetPrecise = [ targetPreciseEph( N_Step+1, 2 ); targetPreciseEph( N_Step+1, 3 ); targetPreciseEph( N_Step+1, 4 ) ]./10^3;
vECIManouverEnd_targetPrecise = [ targetPreciseEph( N_Step+1, 5 ); targetPreciseEph( N_Step+1, 6 ); targetPreciseEph( N_Step+1, 7 ) ]./10^3;


%%%%%%%%%%%%%  HPOP MODEL


%%%%%%%  Experiment 1 Only Thrust effect d = 1s


AuxParam.thrustDuration = 5;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 40;
AuxParam.m       = 40;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.sRad    = 1;
AuxParam.drag    = 1;
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
    
    thisDeviationTime = MCtimeDeviation( experimentIndex );

    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    
    N_Step_Initial = floor((maneuverStartTime + thisDeviationTime) *1/Step); % 
    N_Step_Thrust = round(AuxParam.thrustDuration*1/Step); %  
    N_Step_Final = ceil(maneuverEndTime*1/Step) - N_Step_Initial - N_Step_Thrust; %  
    
    ErronousVelocityChange = AuxParam.velocityChangeLVLH .* MCthrustOutputDeviation( experimentIndex );

    AuxParam.thrustLVLHAcceleration = ErronousVelocityChange./AuxParam.thrustDuration;
    
    RMatThrustDeviation = RotMatXYZEulerDeg(MCthrustDirectionRollDeviationDeg( experimentIndex ), MCthrustDirectionPitchDeviationDeg( experimentIndex ), MCthrustDirectionYawDeviationDeg( experimentIndex ));
    AuxParam.thrustLVLHAcceleration = RMatThrustDeviation * AuxParam.thrustLVLHAcceleration;

    AuxParam.Thrust = 0;
    [initialEph] = ephemeris_v3(Y0, N_Step_Initial, Step);
    currY = [ initialEph(N_Step_Initial+1, 2), initialEph(N_Step_Initial+1, 3), initialEph(N_Step_Initial+1, 4), initialEph(N_Step_Initial+1, 5), initialEph(N_Step_Initial+1, 6), initialEph(N_Step_Initial+1, 7) ];

    AuxParam.Thrust = 1;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime  + MCtimeDeviation( 1 ))/const.DAYSEC);
    [thrustEph] = ephemeris_v3(currY, N_Step_Thrust, Step);
    currY = [ thrustEph(N_Step_Thrust+1, 2), thrustEph(N_Step_Thrust+1, 3), thrustEph(N_Step_Thrust+1, 4), thrustEph(N_Step_Thrust+1, 5), thrustEph(N_Step_Thrust+1, 6), thrustEph(N_Step_Thrust+1, 7) ];

    AuxParam.Thrust = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + AuxParam.thrustDuration + MCtimeDeviation( 1 ))/const.DAYSEC);
    [finalEph] = ephemeris_v3(currY, N_Step_Final, Step);

    Eph = [initialEph; thrustEph; finalEph];

    MC_1_HPOP_PosEnd( experimentIndex, : ) = [Eph(N_Step+3, 2), Eph(N_Step+3, 3), Eph(N_Step+3, 4)]./10^3;
    MC_1_HPOP_VelEnd( experimentIndex, : ) = [Eph(N_Step+3, 5), Eph(N_Step+3, 6), Eph(N_Step+3, 7)]./10^3;

    MC_1_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph(:, 2)'./10^3;
    MC_1_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph(:, 3)'./10^3;
    MC_1_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph(:, 4)'./10^3;
    
end



%% Statistical Analysis

[ QmatECItoLVLH_targetEnd ] = ECIToLVLH( rECIManouverEnd_targetPrecise, vECIManouverEnd_targetPrecise );

absDeviationEndPosHPOP_1 = zeros( MCsampleNum, 1 );
relEndPosECIHPOP_1_chaser = MC_1_HPOP_PosEnd - ones(MCsampleNum,1)*rECIManouverEnd_targetPrecise';
relEndPosLVLHHPOP_1_chaser = zeros(MCsampleNum,3);
for posIter = 1:MCsampleNum
    thisECIPos = relEndPosECIHPOP_1_chaser( posIter, : )';
    thisLVLHPos = QmatECItoLVLH_targetEnd * thisECIPos;
    relEndPosLVLHHPOP_1_chaser( posIter, : ) = thisLVLHPos';
end


relXTrajectoryHPOP_1_chaser = zeros(MCsampleNum,N_Step+1);
relYTrajectoryHPOP_1_chaser = zeros(MCsampleNum,N_Step+1);
relZTrajectoryHPOP_1_chaser = zeros(MCsampleNum,N_Step+1);
for MCIter = 1:MCsampleNum
    thisECIXTrajectory = MC_1_HPOP_ECI_X_Trajectories(MCIter, :);
    thisECIYTrajectory = MC_1_HPOP_ECI_Y_Trajectories(MCIter, :);
    thisECIZTrajectory = MC_1_HPOP_ECI_Z_Trajectories(MCIter, :);
    for trajectoryIter = 1:N_Step+1
        thisECITargetPos = [ targetPreciseEph(trajectoryIter, 2); targetPreciseEph(trajectoryIter, 3); targetPreciseEph(trajectoryIter, 4) ]./10^3;
        thisECITargetVel = [ targetPreciseEph(trajectoryIter, 5); targetPreciseEph(trajectoryIter, 6); targetPreciseEph(trajectoryIter, 7) ]./10^3;
        QmatECItoLVLH_targetThis = ECIToLVLH( thisECITargetPos, thisECITargetVel );
        thisChaserECIPos = [ thisECIXTrajectory( trajectoryIter + 2 ); thisECIYTrajectory( trajectoryIter + 2 ); thisECIZTrajectory( trajectoryIter + 2 ) ];
        thisECIRelPos = thisChaserECIPos - thisECITargetPos;
        thisLVLHRelPos = QmatECItoLVLH_targetThis * thisECIRelPos;
        relXTrajectoryHPOP_1_chaser( MCIter, trajectoryIter ) = thisLVLHRelPos(1);
        relYTrajectoryHPOP_1_chaser( MCIter, trajectoryIter ) = thisLVLHRelPos(2);
        relZTrajectoryHPOP_1_chaser( MCIter, trajectoryIter ) = thisLVLHRelPos(3);
    end
    
end



for dataIndex = 1 : MCsampleNum

    thisEndPosHPOP_1 = MC_1_HPOP_PosEnd( dataIndex, : );

    absDeviationEndPosHPOP_1( dataIndex ) = norm( rECIManouverEnd_targetPrecise - thisEndPosHPOP_1' );

end

absMeanDeviationEndPosHPOP_1 = zeros( MCsampleNum, 1 );

for dataIndex = 1 : MCsampleNum

    absMeanDeviationEndPosHPOP_1( dataIndex ) = mean( absDeviationEndPosHPOP_1( 1 : dataIndex ));

end


%% Plots


figure(1)
hold on
grid on
title('Error in ECI XYZ-Position of MC Simulations')
plot3( relEndPosECIHPOP_1_chaser( :, 1 ), relEndPosECIHPOP_1_chaser( :, 2 ), relEndPosECIHPOP_1_chaser( :, 3 ), '*' )
plot3(0,0,0,'m+', 'linewidth',8)
legend('HPOP 1', 'Goal Position')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
hold off

%%

figure(2)
hold on
grid on
title('Maneuver End Position LVLH')
plot3( relEndPosLVLHHPOP_1_chaser( :, 1 ), relEndPosLVLHHPOP_1_chaser( :, 2 ), relEndPosLVLHHPOP_1_chaser( :, 3 ), '*' )
plot3(0,0,0,'m+', 'linewidth',8)
legend('HPOP 1', 'Goal Position')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
hold off

%%

figure(3)
hold on
grid on
title('Norm of Error in Point of Rendezvous of MC Simulations')
plot( absDeviationEndPosHPOP_1 )
legend('HPOP 1')
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
%[ sx, sy, sz ] = sphere;
%surf( sx*rEarth, sy*rEarth, sz*rEarth, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.05 );
for plotIndex = 1 : MCsampleNum
    plot3( MC_1_HPOP_ECI_X_Trajectories(plotIndex, :), MC_1_HPOP_ECI_Y_Trajectories(plotIndex, :), MC_1_HPOP_ECI_Z_Trajectories(plotIndex, :))
end
plot3( rECIManouverEnd_targetPrecise(1), rECIManouverEnd_targetPrecise(2), rECIManouverEnd_targetPrecise(3), '*k' )
text( rECIManouverEnd_targetPrecise(1), rECIManouverEnd_targetPrecise(2), rECIManouverEnd_targetPrecise(3), 'Chaser Ideal Maneuver End' )
hold off


figure(5)
hold on
grid on
title('Target LVLH Trajectories')
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
for plotIndex = 1 : MCsampleNum
    plot3( relXTrajectoryHPOP_1_chaser(plotIndex, :), relYTrajectoryHPOP_1_chaser(plotIndex, :), relZTrajectoryHPOP_1_chaser(plotIndex, :))
end
plot3(0,0,0,'m+', 'linewidth',8)
hold off

%%

figure(7)
hold on
grid on
title('Mean Norm of End Position Error')
plot( absMeanDeviationEndPosHPOP_1 )
legend('HPOP')
xlabel('Sample')
ylabel('Distance [km]')
hold off


%%

figure(8)
hold on
grid on
histogram(absDeviationEndPosHPOP_1,'Normalization','probability','BinWidth',0.5)
xlabel('Distance [km]')
ylabel('Probability')
title('Normalized Histogram End Position Error')
legend('HPOP')
hold off

%%

figure(9)
hold on
grid on

% Create plots.
t = tiledlayout(3,1);
ax1 = nexttile;

histogram(ax1, relEndPosLVLHHPOP_1_chaser(:,1),'Normalization','probability')
legend(ax1,'X-axis')

ax2 = nexttile;

histogram(ax2, relEndPosLVLHHPOP_1_chaser(:,2),'Normalization','probability')
legend(ax2,'Y-axis')

ax3 = nexttile;

histogram(ax3, relEndPosLVLHHPOP_1_chaser(:,3),'Normalization','probability')
legend(ax3,'Z-axis')


% Link the axes
linkaxes([ax1,ax2, ax3],'x');

title(t,'Normalized Histogram Target LVLH End Position Error')
xlabel(t,'Distance [km]')
ylabel(t,'Probability')

hold off

%%

figure(10)
hold on
grid on
histogram(MCtimeDeviation,'Normalization','probability')
xlabel('Time [s]')
ylabel('Probability')
title('Normalized Histogram Thrust Timing Deviation')
hold off

%%

figure(11)
hold on
grid on

% Create plots.
t = tiledlayout(3,1);
ax1 = nexttile;

histogram(ax1, MCthrustDirectionRollDeviationDeg,'Normalization','probability')
legend(ax1,'Roll')

ax2 = nexttile;

histogram(ax2, MCthrustDirectionPitchDeviationDeg,'Normalization','probability')
legend(ax2,'Pitch')

ax3 = nexttile;

histogram(ax3, MCthrustDirectionYawDeviationDeg,'Normalization','probability')
legend(ax3,'Yaw')


% Link the axes
linkaxes([ax1,ax2, ax3],'x');

title(t,'Normalized Histogram Thrust Direction Deviation')
xlabel(t,'Deviation [Deg]')
ylabel(t,'Probability')

hold off

%%

figure(12)
hold on
grid on
histogram(MCthrustOutputDeviation,'Normalization','probability')
xlabel('Performance Factor')
ylabel('Probability')
title('Normalized Histogram Thrust Output Deviation')
hold off
