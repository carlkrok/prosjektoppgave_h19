%% Clear Workspace
%clc
clear all
format long g
%close all

disp('Workspace Cleared')

%% Load Parameters

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC 

run earthParametersHPOP; 

run SatelliteInitialPositions3;


% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0,...
                  'Thrust',0, 'stepCounter', 0, 'velocityChangeLVLH', zeros(3,1),...
                  'thrustECIAcceleration', zeros(3,1), 'prevTimeStep', 0,...
                  'accelIntegral', zeros(3,1), 'thrustDuration', 0,...
                  'velocityChangeECI', zeros(3,1), 'thrustInitiated', 0,...
                  'thrustLVLHAcceleration', zeros(3,1),...
                  'thrustRotMat', 0, 'thrustECIStartDir', zeros(3,1));

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
Y0_chaser1 = [ r0ECI_chaser1', v0ECI_chaser1'].*10^3;
Y0_chaser2 = [ r0ECI_chaser2', v0ECI_chaser2'].*10^3;
Y0_chaser3 = [ r0ECI_chaser3', v0ECI_chaser3'].*10^3;
Y0_chaser4 = [ r0ECI_chaser4', v0ECI_chaser4'].*10^3;
AuxParam.area_solar = 0.2;
AuxParam.area_drag = 0.1;
AuxParam.mass = 2.0;
AuxParam.Cr = 1.0;
AuxParam.Cd = 2.2;

% epoch
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
%Y0 = ECEF2ECI(Mjd_UTC, Y0); Already in ECI

AuxParam.Mjd_UTC = Mjd_UTC;

disp('Parameters Loaded')

%% Experiment Setup

maneuverEndTime = maneuverStartDelay + maneuverTime;
maneuverStartTime = maneuverStartDelay;
maneuverTimeBuffer = 10;
thrustPrecisionFactor = 10000; 
thrustDuration = 0.001;


Step   = 1;   % [s]
N_Step = round(maneuverEndTime*1/Step); 
N_Step_Initial = round(maneuverStartTime *1/Step);
              
Mjd0   = Mjd_UTC;
AuxParam.Mjd_UTC = Mjd_UTC;

% shorten PC, eopdata, swdata, Cnm, and Snm
num = fix(N_Step*Step/86400)+2;
JD = Mjd_UTC+2400000.5;
pointIter = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
PC = PC(pointIter:pointIter+num,:);
mjd = (floor(Mjd_UTC));
pointIter = find(mjd==eopdata(4,:),1,'first');
eopdata = eopdata(:,pointIter:pointIter+num);
pointIter = find((year==swdata(1,:)) & (mon==swdata(2,:)) & (day==swdata(3,:)),1,'first');
swdata = swdata(:,pointIter-3:pointIter+num);
%Cnm = Cnm(1:AuxParam.n+1,1:AuxParam.n+1);
%Snm = Snm(1:AuxParam.n+1,1:AuxParam.n+1);

orbitType_chaser = orbitType;


%% Monte Carlo Experiment Setup

MCsampleNum = 100;

meanDeviationTimeSetup = 0;
%maxDeviationTimeSetup = 0.1;
stdDeviationTimeSetup = 0.5;%1;

%MCtimeDeviation = meanDeviationTimeSetup - stdDeviationTimeSetup + ( 2 * stdDeviationTimeSetup * rand( MCsampleNum, 1 ) );
%MCtimeDeviation = (meanDeviationTimeSetup-maxDeviationTimeSetup : (2*maxDeviationTimeSetup)/(MCsampleNum-1) : meanDeviationTimeSetup+maxDeviationTimeSetup)';
MC_1_timeDeviation = meanDeviationTimeSetup + stdDeviationTimeSetup .* randn( MCsampleNum, 1 );
MC_2_timeDeviation = meanDeviationTimeSetup + stdDeviationTimeSetup .* randn( MCsampleNum, 1 );
MC_3_timeDeviation = meanDeviationTimeSetup + stdDeviationTimeSetup .* randn( MCsampleNum, 1 );
MC_4_timeDeviation = meanDeviationTimeSetup + stdDeviationTimeSetup .* randn( MCsampleNum, 1 );
%MCtimeDeviation = [0,-2:0.4:2,zeros(1,11)];


meanThrustOutputFactor = 1;
stdThrustOutputUncertaintyFactor = 0.001;

%MCthrustOutputDeviation = meanThrustOutputFactor - stdThrustOutputUncertaintyFactor + ( 2 * stdThrustOutputUncertaintyFactor * rand( MCsampleNum, 1 ) );
%MCthrustOutputDeviation = (meanThrustOutputFactor-stdThrustOutputUncertaintyFactor : (2*stdThrustOutputUncertaintyFactor)/(MCsampleNum-1) : meanThrustOutputFactor+stdThrustOutputUncertaintyFactor)';
MC_1_thrustOutputDeviation = meanThrustOutputFactor + stdThrustOutputUncertaintyFactor...
        .* randn( MCsampleNum, 1 );
MC_2_thrustOutputDeviation = meanThrustOutputFactor + stdThrustOutputUncertaintyFactor...
        .* randn( MCsampleNum, 1 );
MC_3_thrustOutputDeviation = meanThrustOutputFactor + stdThrustOutputUncertaintyFactor...
        .* randn( MCsampleNum, 1 );
MC_4_thrustOutputDeviation = meanThrustOutputFactor + stdThrustOutputUncertaintyFactor...
        .* randn( MCsampleNum, 1 );
%MCthrustOutputDeviation = [1,ones(1,11),0.995:0.001:1.005];


meanThrustDirectionErrorRollDeg = 0;
stdThrustDirectionErrorRollDeg = 0.1;

%MCthrustDirectionRollDeviationDeg = meanThrustDirectionErrorRollDeg - stdThrustDirectionErrorRollDeg + ( 2 * stdThrustDirectionErrorRollDeg * rand( MCsampleNum, 1 ) );
%MCthrustDirectionRollDeviationRad = (meanThrustDirectionErrorRollDeg-maxThrustDirectionErrorRollDeg : (2*maxThrustDirectionErrorRollDeg)/(MCsampleNum-1) : meanThrustDirectionErrorRollDeg+maxThrustDirectionErrorRollDeg)'*pi/180;
MC_1_thrustDirectionRollDeviationDeg = meanThrustDirectionErrorRollDeg + ...
    stdThrustDirectionErrorRollDeg .* randn( MCsampleNum, 1 );
MC_2_thrustDirectionRollDeviationDeg = meanThrustDirectionErrorRollDeg + ...
    stdThrustDirectionErrorRollDeg .* randn( MCsampleNum, 1 );
MC_3_thrustDirectionRollDeviationDeg = meanThrustDirectionErrorRollDeg + ...
    stdThrustDirectionErrorRollDeg .* randn( MCsampleNum, 1 );
MC_4_thrustDirectionRollDeviationDeg = meanThrustDirectionErrorRollDeg + ...
    stdThrustDirectionErrorRollDeg .* randn( MCsampleNum, 1 );
%MCthrustDirectionRollDeviationDeg = [0,-2:0.4:2,zeros(1,22)];

meanThrustDirectionErrorPitchDeg = 0;
stdThrustDirectionErrorPitchDeg = 0.1;

%MCthrustDirectionPitchDeviationDeg = meanThrustDirectionErrorPitchDeg - stdThrustDirectionErrorPitchDeg + ( 2 * stdThrustDirectionErrorPitchDeg * rand( MCsampleNum, 1 ) );
%MCthrustDirectionPitchDeviationRad = (meanThrustDirectionErrorPitchDeg-maxThrustDirectionErrorPitchDeg : (2*maxThrustDirectionErrorPitchDeg)/(MCsampleNum-1) : meanThrustDirectionErrorPitchDeg+maxThrustDirectionErrorPitchDeg)'*pi/180;
MC_1_thrustDirectionPitchDeviationDeg = meanThrustDirectionErrorPitchDeg + stdThrustDirectionErrorPitchDeg ...
    .* randn( MCsampleNum, 1 );
MC_2_thrustDirectionPitchDeviationDeg = meanThrustDirectionErrorPitchDeg + stdThrustDirectionErrorPitchDeg ...
    .* randn( MCsampleNum, 1 );
MC_3_thrustDirectionPitchDeviationDeg = meanThrustDirectionErrorPitchDeg + stdThrustDirectionErrorPitchDeg ...
    .* randn( MCsampleNum, 1 );
MC_4_thrustDirectionPitchDeviationDeg = meanThrustDirectionErrorPitchDeg + stdThrustDirectionErrorPitchDeg ...
    .* randn( MCsampleNum, 1 );
%MCthrustDirectionPitchDeviationDeg = [0,zeros(1,11),-2:0.4:2,zeros(1,11)];

meanThrustDirectionErrorYawDeg = 0;
stdThrustDirectionErrorYawDeg = 0.1;

%MCthrustDirectionYawDeviationDeg = meanThrustDirectionErrorYawDeg - stdThrustDirectionErrorYawDeg + ( 2 * stdThrustDirectionErrorYawDeg * rand( MCsampleNum, 1 ) );
%MCthrustDirectionYawDeviationRad = (meanThrustDirectionErrorYawDeg-maxThrustDirectionErrorYawDeg : (2*maxThrustDirectionErrorYawDeg)/(MCsampleNum-1) : meanThrustDirectionErrorYawDeg+maxThrustDirectionErrorYawDeg)'*pi/180;
MC_1_thrustDirectionYawDeviationDeg = meanThrustDirectionErrorYawDeg + stdThrustDirectionErrorYawDeg...
    .* randn( MCsampleNum, 1 );
MC_2_thrustDirectionYawDeviationDeg = meanThrustDirectionErrorYawDeg + stdThrustDirectionErrorYawDeg...
    .* randn( MCsampleNum, 1 );
MC_3_thrustDirectionYawDeviationDeg = meanThrustDirectionErrorYawDeg + stdThrustDirectionErrorYawDeg...
    .* randn( MCsampleNum, 1 );
MC_4_thrustDirectionYawDeviationDeg = meanThrustDirectionErrorYawDeg + stdThrustDirectionErrorYawDeg...
    .* randn( MCsampleNum, 1 );
%MCthrustDirectionYawDeviationDeg = [0,zeros(1,22),-2:0.4:2];

disp('Experiment Defined')

%% Initial Orbit Determination

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
% AuxParam.n       = 0;
% AuxParam.m       = 0;
% AuxParam.sun     = 0;
% AuxParam.moon    = 0;
% AuxParam.planets = 0;
% AuxParam.sRad    = 0;
% AuxParam.drag    = 0;
% AuxParam.SolidEarthTides = 0;
% AuxParam.OceanTides = 0;
% AuxParam.Relativity = 0;

targetY0 = [ r0ECI_target', v0ECI_target' ].*10^3;

AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
targetInitalEph = ephemeris_v4(targetY0, N_Step_Initial, Step);
targetYManeuverStartPrecise = targetInitalEph( end, 2:end );

% B's ideal position at experiment start
chaser1_Y0 = [ r0ECI_chaser1', v0ECI_chaser1' ].*10^3;
AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
chaser1_InitialEph = ephemeris_v4(chaser1_Y0, N_Step_Initial, Step);
rECIVelocityChangePrecise_chaser1 = chaser1_InitialEph( end, 2:4 )'./10^3;
vECIVelocityChangePrecise_chaser1 = chaser1_InitialEph( end, 5:7 )'./10^3;

AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
chaser1_ThrustStartEph = ephemeris_v4(chaser1_Y0, round((maneuverStartTime - 0.5*thrustDuration) *1/Step), Step);
rECIThrustStartPrecise_chaser1 = chaser1_ThrustStartEph( end, 2:4 )'./10^3;
vECIThrustStartPrecise_chaser1 = chaser1_ThrustStartEph( end, 5:7 )'./10^3;
%
chaser2_Y0 = [ r0ECI_chaser2', v0ECI_chaser2' ].*10^3;
AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
chaser2_InitialEph = ephemeris_v4(chaser2_Y0, N_Step_Initial, Step);
rECIVelocityChangePrecise_chaser2 = chaser2_InitialEph( end, 2:4 )'./10^3;
vECIVelocityChangePrecise_chaser2 = chaser2_InitialEph( end, 5:7 )'./10^3;

AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
chaser2_ThrustStartEph = ephemeris_v4(chaser2_Y0, round((maneuverStartTime - 0.5*thrustDuration) *1/Step), Step);
rECIThrustStartPrecise_chaser2 = chaser2_ThrustStartEph( end, 2:4 )'./10^3;
vECIThrustStartPrecise_chaser2 = chaser2_ThrustStartEph( end, 5:7 )'./10^3;
%
chaser3_Y0 = [ r0ECI_chaser3', v0ECI_chaser3' ].*10^3;
AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
chaser3_InitialEph = ephemeris_v4(chaser3_Y0, N_Step_Initial, Step);
rECIVelocityChangePrecise_chaser3 = chaser3_InitialEph( end, 2:4 )'./10^3;
vECIVelocityChangePrecise_chaser3 = chaser3_InitialEph( end, 5:7 )'./10^3;

AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
chaser3_ThrustStartEph = ephemeris_v4(chaser3_Y0, round((maneuverStartTime - 0.5*thrustDuration) *1/Step), Step);
rECIThrustStartPrecise_chaser3 = chaser3_ThrustStartEph( end, 2:4 )'./10^3;
vECIThrustStartPrecise_chaser3 = chaser3_ThrustStartEph( end, 5:7 )'./10^3;
%
chaser4_Y0 = [ r0ECI_chaser4', v0ECI_chaser4' ].*10^3;
AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
chaser4_InitialEph = ephemeris_v4(chaser4_Y0, N_Step_Initial, Step);
rECIVelocityChangePrecise_chaser4 = chaser4_InitialEph( end, 2:4 )'./10^3;
vECIVelocityChangePrecise_chaser4 = chaser4_InitialEph( end, 5:7 )'./10^3;

AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
chaser4_ThrustStartEph = ephemeris_v4(chaser4_Y0, round((maneuverStartTime - 0.5*thrustDuration) *1/Step), Step);
rECIThrustStartPrecise_chaser4 = chaser4_ThrustStartEph( end, 2:4 )'./10^3;
vECIThrustStartPrecise_chaser4 = chaser4_ThrustStartEph( end, 5:7 )'./10^3;




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


AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
targetManeuverKeplerianEph = ephemeris_v4(targetYManeuverStartPrecise, N_Step - N_Step_Initial, Step);

targetMixedEph = [ targetInitalEph( 1:N_Step_Initial, :); targetManeuverKeplerianEph ];

rECIManouverEnd_targetKeplerian = targetMixedEph( end, 2:4 )'./10^3;
vECIManouverEnd_targetKeplerian = targetMixedEph( end, 5:7 )'./10^3;

QmatECItoLVLH_targetEnd = ECIToLVLH( rECIManouverEnd_targetKeplerian, vECIManouverEnd_targetKeplerian );
QmatLVLHtoECI_targetEnd = QmatECItoLVLH_targetEnd';

chaser1_ECI_endPos = rECIManouverEnd_targetKeplerian + QmatLVLHtoECI_targetEnd*[0;0;chaserTargetDistance];
chaser2_ECI_endPos = rECIManouverEnd_targetKeplerian + QmatLVLHtoECI_targetEnd*[chaserTargetDistance;0;0];
chaser3_ECI_endPos = rECIManouverEnd_targetKeplerian + QmatLVLHtoECI_targetEnd*[0;0;-chaserTargetDistance];
chaser4_ECI_endPos = rECIManouverEnd_targetKeplerian + QmatLVLHtoECI_targetEnd*[-chaserTargetDistance;0;0];



% Required velocity change satellite B
[ deltaVStartECI_chaser1, deltaVEndECI_chaser1, vIntersectOrbit_chaser1 ] = interceptOrbit( rECIVelocityChangePrecise_chaser1, vECIVelocityChangePrecise_chaser1, chaser1_ECI_endPos, vECIManouverEnd_targetKeplerian, maneuverEndTime - maneuverStartTime , orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
[ QmatECItoLVLH_chaserVelocityChange ] = ECIToLVLH( rECIVelocityChangePrecise_chaser1, vECIVelocityChangePrecise_chaser1 );
deltaVStartLVLH_chaser1 = QmatECItoLVLH_chaserVelocityChange * deltaVStartECI_chaser1;
chaser1_velocityChangeECI = deltaVStartECI_chaser1*1000;
chaser1_velocityChangeLVLH = deltaVStartLVLH_chaser1*1000;
[ QmatECItoLVLH_chaserThrusterStart ] = ECIToLVLH( rECIThrustStartPrecise_chaser1, vECIThrustStartPrecise_chaser1 );
QmatLVLHtoECI_chaserThrusterStart = QmatECItoLVLH_chaserThrusterStart';
chaser1_thrustECIStartDir = QmatLVLHtoECI_chaserThrusterStart * (chaser1_velocityChangeLVLH ./ norm(chaser1_velocityChangeLVLH));


[ deltaVStartECI_chaser2, deltaVEndECI_chaser2, vIntersectOrbit_chaser2 ] = interceptOrbit( rECIVelocityChangePrecise_chaser2, vECIVelocityChangePrecise_chaser2, chaser2_ECI_endPos, vECIManouverEnd_targetKeplerian, maneuverEndTime - maneuverStartTime , orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
[ QmatECItoLVLH_chaserVelocityChange ] = ECIToLVLH( rECIVelocityChangePrecise_chaser2, vECIVelocityChangePrecise_chaser2 );
deltaVStartLVLH_chaser2 = QmatECItoLVLH_chaserVelocityChange * deltaVStartECI_chaser2;
chaser2_velocityChangeECI = deltaVStartECI_chaser2*1000;
chaser2_velocityChangeLVLH = deltaVStartLVLH_chaser2*1000;
[ QmatECItoLVLH_chaserThrusterStart ] = ECIToLVLH( rECIThrustStartPrecise_chaser2, vECIThrustStartPrecise_chaser2 );
QmatLVLHtoECI_chaserThrusterStart = QmatECItoLVLH_chaserThrusterStart';
chaser2_thrustECIStartDir = QmatLVLHtoECI_chaserThrusterStart * (chaser2_velocityChangeLVLH ./ norm(chaser2_velocityChangeLVLH));


[ deltaVStartECI_chaser3, deltaVEndECI_chaser3, vIntersectOrbit_chaser3 ] = interceptOrbit( rECIVelocityChangePrecise_chaser3, vECIVelocityChangePrecise_chaser3, chaser3_ECI_endPos, vECIManouverEnd_targetKeplerian, maneuverEndTime - maneuverStartTime , orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
[ QmatECItoLVLH_chaserVelocityChange ] = ECIToLVLH( rECIVelocityChangePrecise_chaser3, vECIVelocityChangePrecise_chaser3 );
deltaVStartLVLH_chaser3 = QmatECItoLVLH_chaserVelocityChange * deltaVStartECI_chaser3;
chaser3_velocityChangeECI = deltaVStartECI_chaser3*1000;
chaser3_velocityChangeLVLH = deltaVStartLVLH_chaser3*1000;
[ QmatECItoLVLH_chaserThrusterStart ] = ECIToLVLH( rECIThrustStartPrecise_chaser3, vECIThrustStartPrecise_chaser3 );
QmatLVLHtoECI_chaserThrusterStart = QmatECItoLVLH_chaserThrusterStart';
chaser3_thrustECIStartDir = QmatLVLHtoECI_chaserThrusterStart * (chaser3_velocityChangeLVLH ./ norm(chaser3_velocityChangeLVLH));


[ deltaVStartECI_chaser4, deltaVEndECI_chaser4, vIntersectOrbit_chaser4 ] = interceptOrbit( rECIVelocityChangePrecise_chaser4, vECIVelocityChangePrecise_chaser4, chaser4_ECI_endPos, vECIManouverEnd_targetKeplerian, maneuverEndTime - maneuverStartTime , orbitType_chaser, muEarth, anomalyErrorTolerance, anomalyMaxIterations );
[ QmatECItoLVLH_chaserVelocityChange ] = ECIToLVLH( rECIVelocityChangePrecise_chaser4, vECIVelocityChangePrecise_chaser4 );
deltaVStartLVLH_chaser4 = QmatECItoLVLH_chaserVelocityChange * deltaVStartECI_chaser4;
chaser4_velocityChangeECI = deltaVStartECI_chaser4*1000;
chaser4_velocityChangeLVLH = deltaVStartLVLH_chaser4*1000;
[ QmatECItoLVLH_chaserThrusterStart ] = ECIToLVLH( rECIThrustStartPrecise_chaser4, vECIThrustStartPrecise_chaser4 );
QmatLVLHtoECI_chaserThrusterStart = QmatECItoLVLH_chaserThrusterStart';
chaser4_thrustECIStartDir = QmatLVLHtoECI_chaserThrusterStart * (chaser4_velocityChangeLVLH ./ norm(chaser4_velocityChangeLVLH));



disp('Initial Orbits Determined')

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
% AuxParam.n       = 0;
% AuxParam.m       = 0;
% AuxParam.sun     = 0;
% AuxParam.moon    = 0;
% AuxParam.planets = 0;
% AuxParam.sRad    = 0;
% AuxParam.drag    = 0;
% AuxParam.SolidEarthTides = 0;
% AuxParam.OceanTides = 0;
% AuxParam.Relativity = 0;

AuxParam.prevTimeStep = 0;
AuxParam.stepCounter = 0;
targetPreciseEph = ephemeris_v4(targetY0, N_Step, Step);
rECIManouverEnd_targetPrecise = targetPreciseEph( end, 2:4 )'./10^3;
vECIManouverEnd_targetPrecise = targetPreciseEph( end, 5:7 )'./10^3;

%%%%%%%%%%%%%  HPOP MODEL


%%%%%%%  Experiment 1 Thrust

AuxParam.Mjd_UTC = Mjd0;
AuxParam.thrustDuration = thrustDuration;
% AuxParam.n       = 40;
% AuxParam.m       = 40;
% AuxParam.sun     = 1;
% AuxParam.moon    = 1;
% AuxParam.planets = 1;
% AuxParam.sRad    = 1;
% AuxParam.drag    = 1;
% AuxParam.SolidEarthTides = 0;
% AuxParam.OceanTides = 0;
% AuxParam.Relativity = 0;


MC_1_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_1_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_1_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_1_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_1_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_1_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_1_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_1_HPOP_ECI_velX = zeros(MCsampleNum, N_Step+1);
MC_1_HPOP_ECI_velY = zeros(MCsampleNum, N_Step+1);
MC_1_HPOP_ECI_velZ = zeros(MCsampleNum, N_Step+1);

MC_2_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_2_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_2_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_2_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_2_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_2_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_2_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_2_HPOP_ECI_velX = zeros(MCsampleNum, N_Step+1);
MC_2_HPOP_ECI_velY = zeros(MCsampleNum, N_Step+1);
MC_2_HPOP_ECI_velZ = zeros(MCsampleNum, N_Step+1);

MC_3_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_3_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_3_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_3_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_3_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_3_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_3_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_3_HPOP_ECI_velX = zeros(MCsampleNum, N_Step+1);
MC_3_HPOP_ECI_velY = zeros(MCsampleNum, N_Step+1);
MC_3_HPOP_ECI_velZ = zeros(MCsampleNum, N_Step+1);

MC_4_HPOP_PosEnd = zeros( MCsampleNum, 3);
MC_4_HPOP_PosStart = zeros( MCsampleNum, 3);
MC_4_HPOP_VelEnd = zeros(MCsampleNum, 3);
MC_4_HPOP_VelStart = zeros(MCsampleNum, 3);
MC_4_HPOP_ECI_X_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_4_HPOP_ECI_Y_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_4_HPOP_ECI_Z_Trajectories = zeros(MCsampleNum, N_Step+1);
MC_4_HPOP_ECI_velX = zeros(MCsampleNum, N_Step+1);
MC_4_HPOP_ECI_velY = zeros(MCsampleNum, N_Step+1);
MC_4_HPOP_ECI_velZ = zeros(MCsampleNum, N_Step+1);


precisionStep = Step / thrustPrecisionFactor;

for experimentIndex = 1 : MCsampleNum
    
    fprintf('Starting MC thrust run %d of %d.\n',experimentIndex,MCsampleNum)

    %%%%%%%%%%%%%  HPOP MODEL
    % propagation
    
    N_Step_Initial = round((maneuverStartTime - 0.5 * maneuverTimeBuffer ) *1/Step); % 
    N_Step_Thrust =  round(maneuverTimeBuffer *1/Step); %  
    N_Step_Final = N_Step - N_Step_Initial - N_Step_Thrust; %  

    
    thisDeviationTime = MC_1_timeDeviation( experimentIndex );
    AuxParam.velocityChangeECI = chaser1_velocityChangeECI;
    AuxParam.velocityChangeLVLH = chaser1_velocityChangeLVLH;
    AuxParam.thrustECIStartDir = chaser1_thrustECIStartDir;
    ErronousVelocityChangeECI = AuxParam.velocityChangeECI .* MC_1_thrustOutputDeviation( experimentIndex );
    AuxParam.thrustECIAcceleration = (ErronousVelocityChangeECI ./ ( AuxParam.thrustDuration  ));
    RMatThrustDeviation = RotMatXYZEulerDeg(MC_1_thrustDirectionRollDeviationDeg( experimentIndex ), MC_1_thrustDirectionPitchDeviationDeg( experimentIndex ), MC_1_thrustDirectionYawDeviationDeg( experimentIndex ));
    AuxParam.thrustECIAcceleration = RMatThrustDeviation * AuxParam.thrustECIAcceleration;
    AuxParam.Thrust = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    AuxParam.Mjd_UTC = Mjd0;
    [initialEph_chaser1] = ephemeris_v4(Y0_chaser1, N_Step_Initial, Step);
    currY_chaser1 = initialEph_chaser1(end, 2:7);
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime  -0.5*maneuverTimeBuffer -0.5*thrustDuration )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_BeforeThrust = round((0.5 * maneuverTimeBuffer - 0.5*thrustDuration + thisDeviationTime) *1/precisionStep);
    [beforeThrustEph_chaser1] = ephemeris_v4(currY_chaser1, N_Step_BeforeThrust, precisionStep);
    currY_chaser1 = [ beforeThrustEph_chaser1(end, 2), beforeThrustEph_chaser1(end, 3), beforeThrustEph_chaser1(end, 4), beforeThrustEph_chaser1(end, 5), beforeThrustEph_chaser1(end, 6), beforeThrustEph_chaser1(end, 7) ];
    AuxParam.Thrust = 1;
    AuxParam.thrustInitiated = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + 0.5*maneuverTimeBuffer -0.5*thrustDuration + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_PrecisionThrust = round((thrustDuration) *1/precisionStep);
    [precisionThrustEph_chaser1] = ephemeris_v4(currY_chaser1, N_Step_PrecisionThrust, precisionStep);
    currY_chaser1 = [ precisionThrustEph_chaser1(end, 2), precisionThrustEph_chaser1(end, 3), precisionThrustEph_chaser1(end, 4), precisionThrustEph_chaser1(end, 5), precisionThrustEph_chaser1(end, 6), precisionThrustEph_chaser1(end, 7) ];
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + 0.5*maneuverTimeBuffer +0.5*thrustDuration + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_AfterThrust = (N_Step_Thrust*thrustPrecisionFactor) - N_Step_PrecisionThrust - N_Step_BeforeThrust; %round((0.5 * maneuverTimeBuffer -0.5*thrustDuration - thisDeviationTime) *1/precisionStep);
    [afterThrustEph_chaser1] = ephemeris_v4(currY_chaser1, N_Step_AfterThrust, precisionStep);
    currY_chaser1 = [ afterThrustEph_chaser1(end, 2), afterThrustEph_chaser1(end, 3), afterThrustEph_chaser1(end, 4), afterThrustEph_chaser1(end, 5), afterThrustEph_chaser1(end, 6), afterThrustEph_chaser1(end, 7) ];
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + maneuverTimeBuffer + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    [finalEph_chaser1] = ephemeris_v4(currY_chaser1, N_Step_Final, Step);
    combinedThrustEph_chaser1 = [beforeThrustEph_chaser1(1:N_Step_BeforeThrust, :); precisionThrustEph_chaser1(1:N_Step_PrecisionThrust, :); afterThrustEph_chaser1(1:N_Step_AfterThrust, :)];
    shortenedThrustEph_chaser1 = combinedThrustEph_chaser1(1:thrustPrecisionFactor:end, :);
    Eph_chaser1 = [initialEph_chaser1(1:N_Step_Initial, :); shortenedThrustEph_chaser1; finalEph_chaser1]; %

    MC_1_HPOP_PosEnd( experimentIndex, : ) = Eph_chaser1( end, 2:4 )./10^3;
    MC_1_HPOP_VelEnd( experimentIndex, : ) = Eph_chaser1( end, 5:7 )./10^3;
    MC_1_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph_chaser1( 1:end , 2)'./10^3;
    MC_1_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph_chaser1( 1:end , 3)'./10^3;
    MC_1_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph_chaser1( 1:end , 4)'./10^3;
    MC_1_HPOP_ECI_velX (experimentIndex, : ) = Eph_chaser1( 1:end , 5 )';
    MC_1_HPOP_ECI_velY (experimentIndex, : ) = Eph_chaser1( 1:end , 6 )';
    MC_1_HPOP_ECI_velZ (experimentIndex, : ) = Eph_chaser1( 1:end , 7 )';
    
    
    thisDeviationTime = MC_2_timeDeviation( experimentIndex );
    AuxParam.velocityChangeECI = chaser2_velocityChangeECI;
    AuxParam.velocityChangeLVLH = chaser2_velocityChangeLVLH;
    AuxParam.thrustECIStartDir = chaser2_thrustECIStartDir;
    ErronousVelocityChangeECI = AuxParam.velocityChangeECI .* MC_1_thrustOutputDeviation( experimentIndex );
    AuxParam.thrustECIAcceleration = (ErronousVelocityChangeECI ./ ( AuxParam.thrustDuration  ));
    RMatThrustDeviation = RotMatXYZEulerDeg(MC_1_thrustDirectionRollDeviationDeg( experimentIndex ), MC_1_thrustDirectionPitchDeviationDeg( experimentIndex ), MC_1_thrustDirectionYawDeviationDeg( experimentIndex ));
    AuxParam.thrustECIAcceleration = RMatThrustDeviation * AuxParam.thrustECIAcceleration;
    AuxParam.Thrust = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    AuxParam.Mjd_UTC = Mjd0;
    [initialEph_chaser2] = ephemeris_v4(Y0_chaser2, N_Step_Initial, Step);
    currY_chaser2 = initialEph_chaser2(end, 2:7);
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime  -0.5*maneuverTimeBuffer -0.5*thrustDuration )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_BeforeThrust = round((0.5 * maneuverTimeBuffer - 0.5*thrustDuration + thisDeviationTime) *1/precisionStep);
    [beforeThrustEph_chaser2] = ephemeris_v4(currY_chaser2, N_Step_BeforeThrust, precisionStep);
    currY_chaser2 = [ beforeThrustEph_chaser2(end, 2), beforeThrustEph_chaser2(end, 3), beforeThrustEph_chaser2(end, 4), beforeThrustEph_chaser2(end, 5), beforeThrustEph_chaser2(end, 6), beforeThrustEph_chaser2(end, 7) ];
    AuxParam.Thrust = 1;
    AuxParam.thrustInitiated = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + 0.5*maneuverTimeBuffer -0.5*thrustDuration + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_PrecisionThrust = round((thrustDuration) *1/precisionStep);
    [precisionThrustEph_chaser2] = ephemeris_v4(currY_chaser2, N_Step_PrecisionThrust, precisionStep);
    currY_chaser2 = [ precisionThrustEph_chaser2(end, 2), precisionThrustEph_chaser2(end, 3), precisionThrustEph_chaser2(end, 4), precisionThrustEph_chaser2(end, 5), precisionThrustEph_chaser2(end, 6), precisionThrustEph_chaser2(end, 7) ];
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + 0.5*maneuverTimeBuffer +0.5*thrustDuration + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_AfterThrust = (N_Step_Thrust*thrustPrecisionFactor) - N_Step_PrecisionThrust - N_Step_BeforeThrust; %round((0.5 * maneuverTimeBuffer -0.5*thrustDuration - thisDeviationTime) *1/precisionStep);
    [afterThrustEph_chaser2] = ephemeris_v4(currY_chaser2, N_Step_AfterThrust, precisionStep);
    currY_chaser2 = [ afterThrustEph_chaser2(end, 2), afterThrustEph_chaser2(end, 3), afterThrustEph_chaser2(end, 4), afterThrustEph_chaser2(end, 5), afterThrustEph_chaser2(end, 6), afterThrustEph_chaser2(end, 7) ];
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + maneuverTimeBuffer + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    [finalEph_chaser2] = ephemeris_v4(currY_chaser2, N_Step_Final, Step);
    combinedThrustEph_chaser2 = [beforeThrustEph_chaser2(1:N_Step_BeforeThrust, :); precisionThrustEph_chaser2(1:N_Step_PrecisionThrust, :); afterThrustEph_chaser2(1:N_Step_AfterThrust, :)];
    shortenedThrustEph_chaser2 = combinedThrustEph_chaser2(1:thrustPrecisionFactor:end, :);
    Eph_chaser2 = [initialEph_chaser2(1:N_Step_Initial, :); shortenedThrustEph_chaser2; finalEph_chaser2]; %

    MC_2_HPOP_PosEnd( experimentIndex, : ) = Eph_chaser2( end, 2:4 )./10^3;
    MC_2_HPOP_VelEnd( experimentIndex, : ) = Eph_chaser2( end, 5:7 )./10^3;
    MC_2_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph_chaser2( 1:end , 2)'./10^3;
    MC_2_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph_chaser2( 1:end , 3)'./10^3;
    MC_2_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph_chaser2( 1:end , 4)'./10^3;
    MC_2_HPOP_ECI_velX (experimentIndex, : ) = Eph_chaser2( 1:end , 5 )';
    MC_2_HPOP_ECI_velY (experimentIndex, : ) = Eph_chaser2( 1:end , 6 )';
    MC_2_HPOP_ECI_velZ (experimentIndex, : ) = Eph_chaser2( 1:end , 7 )';
    
    
    thisDeviationTime = MC_3_timeDeviation( experimentIndex );
    AuxParam.velocityChangeECI = chaser3_velocityChangeECI;
    AuxParam.velocityChangeLVLH = chaser3_velocityChangeLVLH;
    AuxParam.thrustECIStartDir = chaser3_thrustECIStartDir;
    ErronousVelocityChangeECI = AuxParam.velocityChangeECI .* MC_1_thrustOutputDeviation( experimentIndex );
    AuxParam.thrustECIAcceleration = (ErronousVelocityChangeECI ./ ( AuxParam.thrustDuration  ));
    RMatThrustDeviation = RotMatXYZEulerDeg(MC_1_thrustDirectionRollDeviationDeg( experimentIndex ), MC_1_thrustDirectionPitchDeviationDeg( experimentIndex ), MC_1_thrustDirectionYawDeviationDeg( experimentIndex ));
    AuxParam.thrustECIAcceleration = RMatThrustDeviation * AuxParam.thrustECIAcceleration;
    AuxParam.Thrust = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    AuxParam.Mjd_UTC = Mjd0;
    [initialEph_chaser3] = ephemeris_v4(Y0_chaser3, N_Step_Initial, Step);
    currY_chaser3 = initialEph_chaser3(end, 2:7);
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime  -0.5*maneuverTimeBuffer -0.5*thrustDuration )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_BeforeThrust = round((0.5 * maneuverTimeBuffer - 0.5*thrustDuration + thisDeviationTime) *1/precisionStep);
    [beforeThrustEph_chaser3] = ephemeris_v4(currY_chaser3, N_Step_BeforeThrust, precisionStep);
    currY_chaser3 = [ beforeThrustEph_chaser3(end, 2), beforeThrustEph_chaser3(end, 3), beforeThrustEph_chaser3(end, 4), beforeThrustEph_chaser3(end, 5), beforeThrustEph_chaser3(end, 6), beforeThrustEph_chaser3(end, 7) ];
    AuxParam.Thrust = 1;
    AuxParam.thrustInitiated = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + 0.5*maneuverTimeBuffer -0.5*thrustDuration + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_PrecisionThrust = round((thrustDuration) *1/precisionStep);
    [precisionThrustEph_chaser3] = ephemeris_v4(currY_chaser3, N_Step_PrecisionThrust, precisionStep);
    currY_chaser3 = [ precisionThrustEph_chaser3(end, 2), precisionThrustEph_chaser3(end, 3), precisionThrustEph_chaser3(end, 4), precisionThrustEph_chaser3(end, 5), precisionThrustEph_chaser3(end, 6), precisionThrustEph_chaser3(end, 7) ];
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + 0.5*maneuverTimeBuffer +0.5*thrustDuration + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_AfterThrust = (N_Step_Thrust*thrustPrecisionFactor) - N_Step_PrecisionThrust - N_Step_BeforeThrust; %round((0.5 * maneuverTimeBuffer -0.5*thrustDuration - thisDeviationTime) *1/precisionStep);
    [afterThrustEph_chaser3] = ephemeris_v4(currY_chaser3, N_Step_AfterThrust, precisionStep);
    currY_chaser3 = [ afterThrustEph_chaser3(end, 2), afterThrustEph_chaser3(end, 3), afterThrustEph_chaser3(end, 4), afterThrustEph_chaser3(end, 5), afterThrustEph_chaser3(end, 6), afterThrustEph_chaser3(end, 7) ];
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + maneuverTimeBuffer + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    [finalEph_chaser3] = ephemeris_v4(currY_chaser3, N_Step_Final, Step);
    combinedThrustEph_chaser3 = [beforeThrustEph_chaser3(1:N_Step_BeforeThrust, :); precisionThrustEph_chaser3(1:N_Step_PrecisionThrust, :); afterThrustEph_chaser3(1:N_Step_AfterThrust, :)];
    shortenedThrustEph_chaser3 = combinedThrustEph_chaser3(1:thrustPrecisionFactor:end, :);
    Eph_chaser3 = [initialEph_chaser3(1:N_Step_Initial, :); shortenedThrustEph_chaser3; finalEph_chaser3]; %

    MC_3_HPOP_PosEnd( experimentIndex, : ) = Eph_chaser3( end, 2:4 )./10^3;
    MC_3_HPOP_VelEnd( experimentIndex, : ) = Eph_chaser3( end, 5:7 )./10^3;
    MC_3_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph_chaser3( 1:end , 2)'./10^3;
    MC_3_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph_chaser3( 1:end , 3)'./10^3;
    MC_3_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph_chaser3( 1:end , 4)'./10^3;
    MC_3_HPOP_ECI_velX (experimentIndex, : ) = Eph_chaser3( 1:end , 5 )';
    MC_3_HPOP_ECI_velY (experimentIndex, : ) = Eph_chaser3( 1:end , 6 )';
    MC_3_HPOP_ECI_velZ (experimentIndex, : ) = Eph_chaser3( 1:end , 7 )';
    
    
    thisDeviationTime = MC_4_timeDeviation( experimentIndex );
    AuxParam.velocityChangeECI = chaser4_velocityChangeECI;
    AuxParam.velocityChangeLVLH = chaser4_velocityChangeLVLH;
    AuxParam.thrustECIStartDir = chaser4_thrustECIStartDir;
    ErronousVelocityChangeECI = AuxParam.velocityChangeECI .* MC_1_thrustOutputDeviation( experimentIndex );
    AuxParam.thrustECIAcceleration = (ErronousVelocityChangeECI ./ ( AuxParam.thrustDuration  ));
    RMatThrustDeviation = RotMatXYZEulerDeg(MC_1_thrustDirectionRollDeviationDeg( experimentIndex ), MC_1_thrustDirectionPitchDeviationDeg( experimentIndex ), MC_1_thrustDirectionYawDeviationDeg( experimentIndex ));
    AuxParam.thrustECIAcceleration = RMatThrustDeviation * AuxParam.thrustECIAcceleration;
    AuxParam.Thrust = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    AuxParam.Mjd_UTC = Mjd0;
    [initialEph_chaser4] = ephemeris_v4(Y0_chaser4, N_Step_Initial, Step);
    currY_chaser4 = initialEph_chaser4(end, 2:7);
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime  -0.5*maneuverTimeBuffer -0.5*thrustDuration )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_BeforeThrust = round((0.5 * maneuverTimeBuffer - 0.5*thrustDuration + thisDeviationTime) *1/precisionStep);
    [beforeThrustEph_chaser4] = ephemeris_v4(currY_chaser4, N_Step_BeforeThrust, precisionStep);
    currY_chaser4 = [ beforeThrustEph_chaser4(end, 2), beforeThrustEph_chaser4(end, 3), beforeThrustEph_chaser4(end, 4), beforeThrustEph_chaser4(end, 5), beforeThrustEph_chaser4(end, 6), beforeThrustEph_chaser4(end, 7) ];
    AuxParam.Thrust = 1;
    AuxParam.thrustInitiated = 0;
    AuxParam.accelIntegral = zeros(3,1);
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + 0.5*maneuverTimeBuffer -0.5*thrustDuration + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_PrecisionThrust = round((thrustDuration) *1/precisionStep);
    [precisionThrustEph_chaser4] = ephemeris_v4(currY_chaser4, N_Step_PrecisionThrust, precisionStep);
    currY_chaser4 = [ precisionThrustEph_chaser4(end, 2), precisionThrustEph_chaser4(end, 3), precisionThrustEph_chaser4(end, 4), precisionThrustEph_chaser4(end, 5), precisionThrustEph_chaser4(end, 6), precisionThrustEph_chaser4(end, 7) ];
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + 0.5*maneuverTimeBuffer +0.5*thrustDuration + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    N_Step_AfterThrust = (N_Step_Thrust*thrustPrecisionFactor) - N_Step_PrecisionThrust - N_Step_BeforeThrust; %round((0.5 * maneuverTimeBuffer -0.5*thrustDuration - thisDeviationTime) *1/precisionStep);
    [afterThrustEph_chaser4] = ephemeris_v4(currY_chaser4, N_Step_AfterThrust, precisionStep);
    currY_chaser4 = [ afterThrustEph_chaser4(end, 2), afterThrustEph_chaser4(end, 3), afterThrustEph_chaser4(end, 4), afterThrustEph_chaser4(end, 5), afterThrustEph_chaser4(end, 6), afterThrustEph_chaser4(end, 7) ];
    AuxParam.Thrust = 0;
    AuxParam.thrustInitiated = 0;
    AuxParam.Mjd_UTC = Mjd0 + ((maneuverStartTime + maneuverTimeBuffer + thisDeviationTime )/const.DAYSEC);
    AuxParam.prevTimeStep = 0;
    AuxParam.stepCounter = 0;
    [finalEph_chaser4] = ephemeris_v4(currY_chaser4, N_Step_Final, Step);
    combinedThrustEph_chaser4 = [beforeThrustEph_chaser4(1:N_Step_BeforeThrust, :); precisionThrustEph_chaser4(1:N_Step_PrecisionThrust, :); afterThrustEph_chaser4(1:N_Step_AfterThrust, :)];
    shortenedThrustEph_chaser4 = combinedThrustEph_chaser4(1:thrustPrecisionFactor:end, :);
    Eph_chaser4 = [initialEph_chaser4(1:N_Step_Initial, :); shortenedThrustEph_chaser4; finalEph_chaser4]; %

    MC_4_HPOP_PosEnd( experimentIndex, : ) = Eph_chaser4( end, 2:4 )./10^3;
    MC_4_HPOP_VelEnd( experimentIndex, : ) = Eph_chaser4( end, 5:7 )./10^3;
    MC_4_HPOP_ECI_X_Trajectories(experimentIndex, :) = Eph_chaser4( 1:end , 2)'./10^3;
    MC_4_HPOP_ECI_Y_Trajectories(experimentIndex, :) = Eph_chaser4( 1:end , 3)'./10^3;
    MC_4_HPOP_ECI_Z_Trajectories(experimentIndex, :) = Eph_chaser4( 1:end , 4)'./10^3;
    MC_4_HPOP_ECI_velX (experimentIndex, : ) = Eph_chaser4( 1:end , 5 )';
    MC_4_HPOP_ECI_velY (experimentIndex, : ) = Eph_chaser4( 1:end , 6 )';
    MC_4_HPOP_ECI_velZ (experimentIndex, : ) = Eph_chaser4( 1:end , 7 )';
    
    
end




%% Statistical Analysis

disp('Statistical Analysis HPOP1 Started')

[ QmatECItoLVLH_targetEnd ] = ECIToLVLH( rECIManouverEnd_targetPrecise, vECIManouverEnd_targetPrecise );

absDeviationEndPosHPOP_chaser1 = zeros( MCsampleNum, 1 );
absDeviationEndPosLVLH_chaser1 = zeros( MCsampleNum, 1 );
relEndPosECIHPOP_chaser1 = MC_1_HPOP_PosEnd - ones(MCsampleNum,1)*rECIManouverEnd_targetPrecise';
relEndPosLVLHHPOP_chaser1 = zeros(MCsampleNum,3);
for posIter = 1:MCsampleNum
    thisECIPos = relEndPosECIHPOP_chaser1( posIter, : )';
    thisLVLHPos = QmatECItoLVLH_targetEnd * thisECIPos;
    relEndPosLVLHHPOP_chaser1( posIter, : ) = thisLVLHPos';
end
relXTrajectoryHPOP_chaser1 = zeros(MCsampleNum,N_Step+1);
relYTrajectoryHPOP_chaser1 = zeros(MCsampleNum,N_Step+1);
relZTrajectoryHPOP_chaser1 = zeros(MCsampleNum,N_Step+1);
for MCIter = 1:MCsampleNum
    thisECIXTrajectory = MC_1_HPOP_ECI_X_Trajectories(MCIter, :);
    thisECIYTrajectory = MC_1_HPOP_ECI_Y_Trajectories(MCIter, :);
    thisECIZTrajectory = MC_1_HPOP_ECI_Z_Trajectories(MCIter, :);
    for trajectoryIter = 1:N_Step+1
        thisECITargetPos = [ targetPreciseEph(trajectoryIter, 2); targetPreciseEph(trajectoryIter, 3); targetPreciseEph(trajectoryIter, 4) ]./10^3;
        thisECITargetVel = [ targetPreciseEph(trajectoryIter, 5); targetPreciseEph(trajectoryIter, 6); targetPreciseEph(trajectoryIter, 7) ]./10^3;
        QmatECItoLVLH_targetThis = ECIToLVLH( thisECITargetPos, thisECITargetVel );
        thisChaserECIPos = [ thisECIXTrajectory( trajectoryIter ); thisECIYTrajectory( trajectoryIter ); thisECIZTrajectory( trajectoryIter ) ];
        thisECIRelPos = thisChaserECIPos - thisECITargetPos;
        thisLVLHRelPos = QmatECItoLVLH_targetThis * thisECIRelPos;
        relXTrajectoryHPOP_chaser1( MCIter, trajectoryIter ) = thisLVLHRelPos(1);
        relYTrajectoryHPOP_chaser1( MCIter, trajectoryIter ) = thisLVLHRelPos(2);
        relZTrajectoryHPOP_chaser1( MCIter, trajectoryIter ) = thisLVLHRelPos(3);
    end 
end
for dataIndex = 1 : MCsampleNum
    thisEndPosHPOP = MC_1_HPOP_PosEnd( dataIndex, : );
    absDeviationEndPosHPOP_chaser1( dataIndex ) = norm( rECIManouverEnd_targetPrecise - thisEndPosHPOP' );
    absDeviationEndPosLVLH_chaser1( dataIndex ) = norm( relEndPosLVLHHPOP_chaser1( dataIndex, : ) );
end
absMeanDeviationEndPosHPOP_chaser1 = zeros( MCsampleNum, 1 );
for dataIndex = 1 : MCsampleNum
    absMeanDeviationEndPosHPOP_chaser1( dataIndex ) = mean( absDeviationEndPosHPOP_chaser1( 1 : dataIndex ));
end


absDeviationEndPosHPOP_chaser2 = zeros( MCsampleNum, 1 );
absDeviationEndPosLVLH_chaser2 = zeros( MCsampleNum, 1 );
relEndPosECIHPOP_chaser2 = MC_2_HPOP_PosEnd - ones(MCsampleNum,1)*rECIManouverEnd_targetPrecise';
relEndPosLVLHHPOP_chaser2 = zeros(MCsampleNum,3);
for posIter = 1:MCsampleNum
    thisECIPos = relEndPosECIHPOP_chaser2( posIter, : )';
    thisLVLHPos = QmatECItoLVLH_targetEnd * thisECIPos;
    relEndPosLVLHHPOP_chaser2( posIter, : ) = thisLVLHPos';
end
relXTrajectoryHPOP_chaser2 = zeros(MCsampleNum,N_Step+1);
relYTrajectoryHPOP_chaser2 = zeros(MCsampleNum,N_Step+1);
relZTrajectoryHPOP_chaser2 = zeros(MCsampleNum,N_Step+1);
for MCIter = 1:MCsampleNum
    thisECIXTrajectory = MC_2_HPOP_ECI_X_Trajectories(MCIter, :);
    thisECIYTrajectory = MC_2_HPOP_ECI_Y_Trajectories(MCIter, :);
    thisECIZTrajectory = MC_2_HPOP_ECI_Z_Trajectories(MCIter, :);
    for trajectoryIter = 1:N_Step+1
        thisECITargetPos = [ targetPreciseEph(trajectoryIter, 2); targetPreciseEph(trajectoryIter, 3); targetPreciseEph(trajectoryIter, 4) ]./10^3;
        thisECITargetVel = [ targetPreciseEph(trajectoryIter, 5); targetPreciseEph(trajectoryIter, 6); targetPreciseEph(trajectoryIter, 7) ]./10^3;
        QmatECItoLVLH_targetThis = ECIToLVLH( thisECITargetPos, thisECITargetVel );
        thisChaserECIPos = [ thisECIXTrajectory( trajectoryIter ); thisECIYTrajectory( trajectoryIter ); thisECIZTrajectory( trajectoryIter ) ];
        thisECIRelPos = thisChaserECIPos - thisECITargetPos;
        thisLVLHRelPos = QmatECItoLVLH_targetThis * thisECIRelPos;
        relXTrajectoryHPOP_chaser2( MCIter, trajectoryIter ) = thisLVLHRelPos(1);
        relYTrajectoryHPOP_chaser2( MCIter, trajectoryIter ) = thisLVLHRelPos(2);
        relZTrajectoryHPOP_chaser2( MCIter, trajectoryIter ) = thisLVLHRelPos(3);
    end 
end
for dataIndex = 1 : MCsampleNum
    thisEndPosHPOP = MC_2_HPOP_PosEnd( dataIndex, : );
    absDeviationEndPosHPOP_chaser2( dataIndex ) = norm( rECIManouverEnd_targetPrecise - thisEndPosHPOP' );
    absDeviationEndPosLVLH_chaser2( dataIndex ) = norm( relEndPosLVLHHPOP_chaser2( dataIndex, : ) );
end
absMeanDeviationEndPosHPOP_chaser2 = zeros( MCsampleNum, 1 );
for dataIndex = 1 : MCsampleNum
    absMeanDeviationEndPosHPOP_chaser2( dataIndex ) = mean( absDeviationEndPosHPOP_chaser2( 1 : dataIndex ));
end


absDeviationEndPosHPOP_chaser3 = zeros( MCsampleNum, 1 );
absDeviationEndPosLVLH_chaser3 = zeros( MCsampleNum, 1 );
relEndPosECIHPOP_chaser3 = MC_3_HPOP_PosEnd - ones(MCsampleNum,1)*rECIManouverEnd_targetPrecise';
relEndPosLVLHHPOP_chaser3 = zeros(MCsampleNum,3);
for posIter = 1:MCsampleNum
    thisECIPos = relEndPosECIHPOP_chaser3( posIter, : )';
    thisLVLHPos = QmatECItoLVLH_targetEnd * thisECIPos;
    relEndPosLVLHHPOP_chaser3( posIter, : ) = thisLVLHPos';
end
relXTrajectoryHPOP_chaser3 = zeros(MCsampleNum,N_Step+1);
relYTrajectoryHPOP_chaser3 = zeros(MCsampleNum,N_Step+1);
relZTrajectoryHPOP_chaser3 = zeros(MCsampleNum,N_Step+1);
for MCIter = 1:MCsampleNum
    thisECIXTrajectory = MC_3_HPOP_ECI_X_Trajectories(MCIter, :);
    thisECIYTrajectory = MC_3_HPOP_ECI_Y_Trajectories(MCIter, :);
    thisECIZTrajectory = MC_3_HPOP_ECI_Z_Trajectories(MCIter, :);
    for trajectoryIter = 1:N_Step+1
        thisECITargetPos = [ targetPreciseEph(trajectoryIter, 2); targetPreciseEph(trajectoryIter, 3); targetPreciseEph(trajectoryIter, 4) ]./10^3;
        thisECITargetVel = [ targetPreciseEph(trajectoryIter, 5); targetPreciseEph(trajectoryIter, 6); targetPreciseEph(trajectoryIter, 7) ]./10^3;
        QmatECItoLVLH_targetThis = ECIToLVLH( thisECITargetPos, thisECITargetVel );
        thisChaserECIPos = [ thisECIXTrajectory( trajectoryIter ); thisECIYTrajectory( trajectoryIter ); thisECIZTrajectory( trajectoryIter ) ];
        thisECIRelPos = thisChaserECIPos - thisECITargetPos;
        thisLVLHRelPos = QmatECItoLVLH_targetThis * thisECIRelPos;
        relXTrajectoryHPOP_chaser3( MCIter, trajectoryIter ) = thisLVLHRelPos(1);
        relYTrajectoryHPOP_chaser3( MCIter, trajectoryIter ) = thisLVLHRelPos(2);
        relZTrajectoryHPOP_chaser3( MCIter, trajectoryIter ) = thisLVLHRelPos(3);
    end 
end
for dataIndex = 1 : MCsampleNum
    thisEndPosHPOP = MC_3_HPOP_PosEnd( dataIndex, : );
    absDeviationEndPosHPOP_chaser3( dataIndex ) = norm( rECIManouverEnd_targetPrecise - thisEndPosHPOP' );
    absDeviationEndPosLVLH_chaser3( dataIndex ) = norm( relEndPosLVLHHPOP_chaser3( dataIndex, : ) );
end
absMeanDeviationEndPosHPOP_chaser3 = zeros( MCsampleNum, 1 );
for dataIndex = 1 : MCsampleNum
    absMeanDeviationEndPosHPOP_chaser3( dataIndex ) = mean( absDeviationEndPosHPOP_chaser3( 1 : dataIndex ));
end


absDeviationEndPosHPOP_chaser4 = zeros( MCsampleNum, 1 );
absDeviationEndPosLVLH_chaser4 = zeros( MCsampleNum, 1 );
relEndPosECIHPOP_chaser4 = MC_4_HPOP_PosEnd - ones(MCsampleNum,1)*rECIManouverEnd_targetPrecise';
relEndPosLVLHHPOP_chaser4 = zeros(MCsampleNum,3);
for posIter = 1:MCsampleNum
    thisECIPos = relEndPosECIHPOP_chaser4( posIter, : )';
    thisLVLHPos = QmatECItoLVLH_targetEnd * thisECIPos;
    relEndPosLVLHHPOP_chaser4( posIter, : ) = thisLVLHPos';
end
relXTrajectoryHPOP_chaser4 = zeros(MCsampleNum,N_Step+1);
relYTrajectoryHPOP_chaser4 = zeros(MCsampleNum,N_Step+1);
relZTrajectoryHPOP_chaser4 = zeros(MCsampleNum,N_Step+1);
for MCIter = 1:MCsampleNum
    thisECIXTrajectory = MC_4_HPOP_ECI_X_Trajectories(MCIter, :);
    thisECIYTrajectory = MC_4_HPOP_ECI_Y_Trajectories(MCIter, :);
    thisECIZTrajectory = MC_4_HPOP_ECI_Z_Trajectories(MCIter, :);
    for trajectoryIter = 1:N_Step+1
        thisECITargetPos = [ targetPreciseEph(trajectoryIter, 2); targetPreciseEph(trajectoryIter, 3); targetPreciseEph(trajectoryIter, 4) ]./10^3;
        thisECITargetVel = [ targetPreciseEph(trajectoryIter, 5); targetPreciseEph(trajectoryIter, 6); targetPreciseEph(trajectoryIter, 7) ]./10^3;
        QmatECItoLVLH_targetThis = ECIToLVLH( thisECITargetPos, thisECITargetVel );
        thisChaserECIPos = [ thisECIXTrajectory( trajectoryIter ); thisECIYTrajectory( trajectoryIter ); thisECIZTrajectory( trajectoryIter ) ];
        thisECIRelPos = thisChaserECIPos - thisECITargetPos;
        thisLVLHRelPos = QmatECItoLVLH_targetThis * thisECIRelPos;
        relXTrajectoryHPOP_chaser4( MCIter, trajectoryIter ) = thisLVLHRelPos(1);
        relYTrajectoryHPOP_chaser4( MCIter, trajectoryIter ) = thisLVLHRelPos(2);
        relZTrajectoryHPOP_chaser4( MCIter, trajectoryIter ) = thisLVLHRelPos(3);
    end 
end
for dataIndex = 1 : MCsampleNum
    thisEndPosHPOP = MC_4_HPOP_PosEnd( dataIndex, : );
    absDeviationEndPosHPOP_chaser4( dataIndex ) = norm( rECIManouverEnd_targetPrecise - thisEndPosHPOP' );
    absDeviationEndPosLVLH_chaser4( dataIndex ) = norm( relEndPosLVLHHPOP_chaser4( dataIndex, : ) );
end
absMeanDeviationEndPosHPOP_chaser4 = zeros( MCsampleNum, 1 );
for dataIndex = 1 : MCsampleNum
    absMeanDeviationEndPosHPOP_chaser4( dataIndex ) = mean( absDeviationEndPosHPOP_chaser4( 1 : dataIndex ));
end


%% Plots



figure(22);
set(gca,'FontSize',30)
set(gcf,'renderer','Painters','Position', [10 10 1500 900])
hold on
grid on
%axis equal
title('Maneuver End Position LVLH')
plot3(0,0,0,'c+', 'linewidth',8)
plot3( relEndPosLVLHHPOP_chaser1( :, 1 ).*10^3, relEndPosLVLHHPOP_chaser1( :, 2 ).*10^3, relEndPosLVLHHPOP_chaser1( :, 3 ).*10^3, '*' )
plot3( relEndPosLVLHHPOP_chaser2( :, 1 ).*10^3, relEndPosLVLHHPOP_chaser2( :, 2 ).*10^3, relEndPosLVLHHPOP_chaser2( :, 3 ).*10^3, '*' )
plot3( relEndPosLVLHHPOP_chaser3( :, 1 ).*10^3, relEndPosLVLHHPOP_chaser3( :, 2 ).*10^3, relEndPosLVLHHPOP_chaser3( :, 3 ).*10^3, '*' )
plot3( relEndPosLVLHHPOP_chaser4( :, 1 ).*10^3, relEndPosLVLHHPOP_chaser4( :, 2 ).*10^3, relEndPosLVLHHPOP_chaser4( :, 3 ).*10^3, '*' )
legend('Target Satellite', 'Chaser 1', 'Chaser 2', 'Chaser 3', 'Chaser 4','Location','northeast')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
view([90 0])
% ylim([-23500 -9500])
% xlim([0 1800])
% zlim([4000 5400])
hold off

%%

figure(5)
hold on
set(gca,'FontSize',30)
%set(gcf,'renderer','Painters','Position', [10 10 1200 600])
grid on
title('Formation maneuver in target LVLH frame')
axis equal
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
h5 = plot3( relEndPosLVLHHPOP_chaser1( :, 1 ).*10^3, relEndPosLVLHHPOP_chaser1( :, 2 ).*10^3, relEndPosLVLHHPOP_chaser1( :, 3 ).*10^3, '*', 'Color', '#0072BD' );
h6 = plot3( relEndPosLVLHHPOP_chaser2( :, 1 ).*10^3, relEndPosLVLHHPOP_chaser2( :, 2 ).*10^3, relEndPosLVLHHPOP_chaser2( :, 3 ).*10^3, '*', 'Color', '#D95319' );
h7 = plot3( relEndPosLVLHHPOP_chaser3( :, 1 ).*10^3, relEndPosLVLHHPOP_chaser3( :, 2 ).*10^3, relEndPosLVLHHPOP_chaser3( :, 3 ).*10^3, '*', 'Color', '#EDB120' );
h8 = plot3( relEndPosLVLHHPOP_chaser4( :, 1 ).*10^3, relEndPosLVLHHPOP_chaser4( :, 2 ).*10^3, relEndPosLVLHHPOP_chaser4( :, 3 ).*10^3, '*', 'Color', '#77AC30' );
h0 = plot3(0,0,0,'r+', 'linewidth',8)
for plotIndex = 1 : MCsampleNum
    h1 = plot3( relXTrajectoryHPOP_chaser1(plotIndex, end-120:end).*10^3, relYTrajectoryHPOP_chaser1(plotIndex, end-120:end).*10^3, relZTrajectoryHPOP_chaser1(plotIndex, end-120:end).*10^3, 'Color', '#0072BD');
    h1.Color(4) = 0.2;
end
for plotIndex = 1 : MCsampleNum
    h2 = plot3( relXTrajectoryHPOP_chaser2(plotIndex, end-120:end).*10^3, relYTrajectoryHPOP_chaser2(plotIndex, end-120:end).*10^3, relZTrajectoryHPOP_chaser2(plotIndex, end-120:end).*10^3, 'Color', '#D95319');
    h2.Color(4) = 0.2;
end
for plotIndex = 1 : MCsampleNum
    h3 = plot3( relXTrajectoryHPOP_chaser3(plotIndex, end-120:end).*10^3, relYTrajectoryHPOP_chaser3(plotIndex, end-120:end).*10^3, relZTrajectoryHPOP_chaser3(plotIndex, end-120:end).*10^3, 'Color', '#EDB120');
    h3.Color(4) = 0.2;
end
for plotIndex = 1 : MCsampleNum
    h4 = plot3( relXTrajectoryHPOP_chaser4(plotIndex, end-120:end).*10^3, relYTrajectoryHPOP_chaser4(plotIndex, end-120:end).*10^3, relZTrajectoryHPOP_chaser4(plotIndex, end-120:end).*10^3, 'Color', '#77AC30');
    h4.Color(4) = 0.2;
end
%view([30 5])
% ylim([-0.05 0.05])
% xlim([-0.05 0.05])
% zlim([-0.05 0.05])
legend([h0,h1,h2,h3,h4],'Target Satellite', 'Chaser 1', 'Chaser 2', 'Chaser 3', 'Chaser 4','Location','eastoutside')
hold off


%%


% 
% xChaserDeviation = zeros(length(relXTrajectoryHPOP_1_chaser),MCsampleNum);
% for pointIter = 1:length(relXTrajectoryHPOP_1_chaser)
%     thisZeroPoint = relXTrajectoryHPOP_1_chaser(1,pointIter);
%     for chaserIter = 2:MCsampleNum
%         xChaserDeviation(pointIter,chaserIter) = relXTrajectoryHPOP_1_chaser(chaserIter,pointIter) - thisZeroPoint;
%     end
% end
% 
% yChaserDeviation = zeros(length(relYTrajectoryHPOP_1_chaser),MCsampleNum);
% for pointIter = 1:length(relYTrajectoryHPOP_1_chaser)
%     thisZeroPoint = relYTrajectoryHPOP_1_chaser(1,pointIter);
%     for chaserIter = 2:MCsampleNum
%         yChaserDeviation(pointIter,chaserIter) = relYTrajectoryHPOP_1_chaser(chaserIter,pointIter) - thisZeroPoint;
%     end
% end
% 
% zChaserDeviation = zeros(length(relZTrajectoryHPOP_1_chaser),MCsampleNum);
% for pointIter = 1:length(relZTrajectoryHPOP_1_chaser)
%     thisZeroPoint = relZTrajectoryHPOP_1_chaser(1,pointIter);
%     for chaserIter = 2:MCsampleNum
%         zChaserDeviation(pointIter,chaserIter) = relZTrajectoryHPOP_1_chaser(chaserIter,pointIter) - thisZeroPoint;
%     end
% end


%%

% figure(71)
% set(gca,'FontSize',30);
%set(gcf,'renderer','Painters')
% hold on
% grid on
% title('Chaser X-axis Deviation Target LVLH')
% xlabel('Time [s]')
% ylabel('Deviation [km]')
% for plotIndex = 2 : 12
%     h1 = plot( 1:length(relXTrajectoryHPOP_1_chaser), xChaserDeviation(:, plotIndex),'r','LineWidth',3);
% end
% for plotIndex = 13 : 23
%     h2 = plot( 1:length(relXTrajectoryHPOP_1_chaser), xChaserDeviation(:, plotIndex),'g','LineWidth',3);
% end
% for plotIndex = 24 : 34
%     h3 = plot( 1:length(relXTrajectoryHPOP_1_chaser), xChaserDeviation(:, plotIndex),'b','LineWidth',3);
% end
% legend([h1,h2,h3],'Roll','Pitch','Yaw','Location','northwest')
% hold off
% 
% figure(72)
% set(gca,'FontSize',30);
%set(gcf,'renderer','Painters')
% hold on
% grid on
% title('Chaser Y-axis Deviation Target LVLH')
% xlabel('Time [s]')
% ylabel('Deviation [km]')
% for plotIndex = 2 : 12
%     h1 = plot( 1:length(relYTrajectoryHPOP_1_chaser), yChaserDeviation(:, plotIndex),'r','LineWidth',3);
% end
% for plotIndex = 13 : 23
%     h2 = plot( 1:length(relYTrajectoryHPOP_1_chaser), yChaserDeviation(:, plotIndex),'g','LineWidth',3);
% end
% for plotIndex = 24 : 34
%     h3 = plot( 1:length(relYTrajectoryHPOP_1_chaser), yChaserDeviation(:, plotIndex),'b','LineWidth',3);
% end
% legend([h1,h2,h3],'Roll','Pitch','Yaw','Location','northwest')
% hold off
% 
% figure(73)
% set(gca,'FontSize',30);
%set(gcf,'renderer','Painters')
% hold on
% grid on
% title('Chaser Z-axis Deviation Target LVLH')
% xlabel('Time [s]')
% ylabel('Deviation [km]')
% for plotIndex = 2 : 12
%     h1 = plot( 1:length(relZTrajectoryHPOP_1_chaser), zChaserDeviation(:, plotIndex),'r','LineWidth',3);
% end
% for plotIndex = 13 : 23
%     h2 = plot( 1:length(relZTrajectoryHPOP_1_chaser), zChaserDeviation(:, plotIndex),'g','LineWidth',3);
% end
% for plotIndex = 24 : 34
%     h3 = plot( 1:length(relZTrajectoryHPOP_1_chaser), zChaserDeviation(:, plotIndex),'b','LineWidth',3);
% end
% legend([h1,h2,h3],'Roll','Pitch','Yaw','Location','northwest')
% hold off

%%

% 
% figure(71)
% set(gca,'FontSize',30);
% set(gcf,'renderer','Painters')
% hold on
% grid on
% title('Chaser X-axis Deviation Target LVLH')
% xlabel('Time [s]')
% ylabel('Deviation [km]')
% for plotIndex = 2 : 12
%     h1 = plot( 1:length(relXTrajectoryHPOP_1_chaser), xChaserDeviation(:, plotIndex),'r','LineWidth',3);
% end
% for plotIndex = 13 : 23
%     h2 = plot( 1:length(relXTrajectoryHPOP_1_chaser), xChaserDeviation(:, plotIndex),'b','LineWidth',3);
% end
% legend([h1,h2],'Timing','Output')
% hold off
% 
% figure(72)
% set(gca,'FontSize',30);
% set(gcf,'renderer','Painters')
% hold on
% grid on
% title('Chaser Y-axis Deviation Target LVLH')
% xlabel('Time [s]')
% ylabel('Deviation [km]')
% for plotIndex = 2 : 12
%     h1 = plot( 1:length(relYTrajectoryHPOP_1_chaser), yChaserDeviation(:, plotIndex),'r','LineWidth',3);
% end
% for plotIndex = 13 : 23
%     h2 = plot( 1:length(relYTrajectoryHPOP_1_chaser), yChaserDeviation(:, plotIndex),'b','LineWidth',3);
% end
% legend([h1,h2],'Timing','Output')
% hold off
% 
% figure(73)
% set(gca,'FontSize',30);
% set(gcf,'renderer','Painters')
% hold on
% grid on
% title('Chaser Z-axis Deviation Target LVLH')
% xlabel('Time [s]')
% ylabel('Deviation [km]')
% for plotIndex = 2 : 12
%     h1 = plot( 1:length(relZTrajectoryHPOP_1_chaser), zChaserDeviation(:, plotIndex),'r','LineWidth',3);
% end
% for plotIndex = 13 : 23
%     h2 = plot( 1:length(relZTrajectoryHPOP_1_chaser), zChaserDeviation(:, plotIndex),'b','LineWidth',3);
% end
% legend([h1,h2],'Timing','Output')
% hold off
% 

