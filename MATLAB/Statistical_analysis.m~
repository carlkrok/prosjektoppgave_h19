
x1 = relEndPosLVLHHPOP_1_chaser.*10^3;
x2 = absDeviationEndPosLVLH_1.*10^3;


%%

mean1 = mean(x1)
mean2 = mean(x2)

%%

var1 = var(x1)
var2 = var(x2)

%%

skew1 = skewness(x1)
skew2 = skewness(x2)

%%

kurt1 = kurtosis(x1)
kurt2 = kurtosis(x2)

%% 

cov1 = cov(x1)


%%

estimated_dist = mvnrnd(mean1,cov1,length(x1));

estimatedPDFVec = mvnpdf(x1,mean1,cov1);


%%

%nBins = 150;
%figure(1)
%hist1 = histogram(x1, nBins, 'Normalization', 'pdf')
%edges1 = hist1.BinEdges;
%figure(2)
%histEst = histogram(estimated_dist, 'BinEdges', edges1, 'Normalization', 'pdf')
%vect1 = hist1.BinCounts;
%vectEst = histEst.BinCounts;


binSizeCube = 100; % Corresponds to Value x Value x Value meters(^3)

maxXSampleVal = max(x1(:,1));
maxYSampleVal = max(x1(:,2));
maxZSampleVal = max(x1(:,3));

minXSampleVal = min(x1(:,1));
minYSampleVal = min(x1(:,2));
minZSampleVal = min(x1(:,3));

xStart = floor(minXSampleVal/binSizeCube)*binSizeCube-0.5*binSizeCube;
xEnd = ceil(maxXSampleVal/binSizeCube)*binSizeCube+0.5*binSizeCube;

yStart = floor(minYSampleVal/binSizeCube)*binSizeCube-0.5*binSizeCube;
yEnd = ceil(maxYSampleVal/binSizeCube)*binSizeCube+0.5*binSizeCube;

zStart = floor(minZSampleVal/binSizeCube)*binSizeCube-0.5*binSizeCube;
zEnd = ceil(maxZSampleVal/binSizeCube)*binSizeCube+0.5*binSizeCube;


numXBins = (xEnd - xStart) / binSizeCube;
numYBins = (yEnd - yStart) / binSizeCube;
numZBins = (zEnd - zStart) / binSizeCube;


binVector = zeros(numXBins, numYBins, numZBins);


for thisPoint = x1'
    
    thisXIndex = ceil((thisPoint(1) - xStart)/binSizeCube);
    thisYIndex = ceil((thisPoint(2) - yStart)/binSizeCube);
    thisZIndex = ceil((thisPoint(3) - zStart)/binSizeCube);
    
    binVector(thisXIndex, thisYIndex, thisZIndex) = binVector(thisXIndex, thisYIndex, thisZIndex) + 1;
    
end

binVectorNormalized = binVector./length(x1);


%%


empiricalProbVec = zeros(length(x1),1);
estimatedProbVec = zeros(length(x1),1);

for pointIndex = 1:length(x1)
    
    thisPoint = x1(pointIndex,:);
    
    thisXIndex = ceil((thisPoint(1) - xStart)/binSizeCube);
    thisYIndex = ceil((thisPoint(2) - yStart)/binSizeCube);
    thisZIndex = ceil((thisPoint(3) - zStart)/binSizeCube);
    
    xLow = xStart + (thisXIndex - 1)*binSizeCube;
    xHigh = xLow + binSizeCube;
    
    yLow = yStart + (thisYIndex - 1)*binSizeCube;
    yHigh = yLow + binSizeCube;
    
    zLow = zStart + (thisZIndex - 1)*binSizeCube;
    zHigh = zLow + binSizeCube;
    
    lb = [xLow; yLow; zLow];
    ub = [xHigh; yHigh; zHigh];
    
    empiricalProbVec(pointIndex) = binVectorNormalized(thisXIndex, thisYIndex, thisZIndex);
    estimatedProbVec(pointIndex) = mvncdf(lb,ub,mean1',cov1);
end


%%

abs(sum(sum(sum(binVectorNormalized))) - 1)


%%

KL_div = sum(empiricalProbVec .* (log2(empiricalProbVec)-log2(estimatedProbVec)))

%%


[coeff,score,latent,tsquared,explained] = pca(x1);

coeff
explained

%%

coeffLarge = coeff*800;

plotHandles = zeros(6,1);

figure(80)
set(gca,'FontSize',20)
set(gcf,'renderer','Painters','Position', [10 10 1000 700])
hold on
axis equal
grid on
title({'Experiment 2 Part 1:','Chaser final position in target LVLH frame'})

a = scatter3(x1(:,1), x1(:,2), x1(:,3), 1,'*');
a.MarkerEdgeAlpha = 0.2;
a.MarkerFaceAlpha = 0.2;
plotHandles(1) = plot3([mean1(1), mean1(1)+coeffLarge(1,1)], [mean1(2), mean1(2)+coeffLarge(2,1)], [mean1(3), mean1(3)+coeffLarge(3,1)],'LineWidth',2);
plotHandles(2) = plot3([mean1(1), mean1(1)+coeffLarge(1,2)], [mean1(2), mean1(2)+coeffLarge(2,2)], [mean1(3), mean1(3)+coeffLarge(3,2)],'LineWidth',2);
plotHandles(3) = plot3([mean1(1), mean1(1)+coeffLarge(1,3)], [mean1(2), mean1(2)+coeffLarge(2,3)], [mean1(3), mean1(3)+coeffLarge(3,3)],'LineWidth',2);
%plotHandles(4) = plot_gaussian_ellipsoid(mean1, cov1, 3);
%b = plot_gaussian_ellipsoid(mean1(1:2), cov1(1:2,1:2), 3);
b = plot_gaussian_2d_xy(mean1, cov1, 3);
c = plot_gaussian_2d_xz(mean1, cov1, 3);
d = plot_gaussian_2d_yz(mean1, cov1, 3);
legend([plotHandles(1:3); b; c; d],'PC1','PC2','PC3','3\sigma_{xy}', '3\sigma_{xz}', '3\sigma_{yz}','Location','eastoutside')
%legend([plotHandles(1:3); b],'PC1','PC2','PC3','3\sigma_{xy}')
%legend([plotHandles(1:3); c],'PC1','PC2','PC3','3\sigma_{xz}')
%legend([plotHandles(1:3); d],'PC1','PC2','PC3','3\sigma_{yz}')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
%a.MarkerFaceAlpha = 0.2;

view([-45 13])
%ylim([9000 14000])
%xlim([0 500])
%zlim([-2100 -1650])
ylim([-6000 6000])
xlim([-1000 1000])
zlim([-600 600])

hold off


%%


chaserTrajNormDev = zeros(MCsampleNum, length(MC_1_HPOP_ECI_X_Trajectories));

for timeIter = 1:length(MC_1_HPOP_ECI_X_Trajectories)
    
    thisPointVector = [MC_1_HPOP_ECI_X_Trajectories(:,timeIter),...
        MC_1_HPOP_ECI_Y_Trajectories(:,timeIter),...
        MC_1_HPOP_ECI_Z_Trajectories(:,timeIter)];
    
    thisAvrgPoint = mean(thisPointVector);
    
    for chaserIter = 1:MCsampleNum
        chaserTrajNormDev(chaserIter,timeIter) = norm(thisPointVector(chaserIter,:) - thisAvrgPoint);
    end
end


figure(91)
hold on
grid on
set(gca,'FontSize',20)
set(gcf,'renderer','Painters')
title('Deviations from Chaser Average Trajectory')
for chaserIter = 1:MCsampleNum
    plot(1:length(MC_1_HPOP_ECI_X_Trajectories),chaserTrajNormDev(chaserIter,:));
end
xlabel('Time [s]')
ylabel('Deviation [km]')
hold off


%%


chaserTrajNormDevX = zeros(MCsampleNum, length(MC_1_HPOP_ECI_X_Trajectories));
chaserTrajNormDevY = zeros(MCsampleNum, length(MC_1_HPOP_ECI_Y_Trajectories));
chaserTrajNormDevZ = zeros(MCsampleNum, length(MC_1_HPOP_ECI_Z_Trajectories));

for timeIter = 1:length(MC_1_HPOP_ECI_X_Trajectories)
    
    thisPointVector = [MC_1_HPOP_ECI_X_Trajectories(:,timeIter),...
        MC_1_HPOP_ECI_Y_Trajectories(:,timeIter),...
        MC_1_HPOP_ECI_Z_Trajectories(:,timeIter)];
    
    thisAvrgPoint = mean(thisPointVector);
    
    thisVelVector = [MC_1_HPOP_ECI_velX(:,timeIter),...
        MC_1_HPOP_ECI_velY(:,timeIter),...
        MC_1_HPOP_ECI_velZ(:,timeIter)];
    
    thisAvrgVel = mean(thisVelVector).*10^-3;
    
    QmatECItoLVLH = ECIToLVLH( thisAvrgPoint', thisAvrgVel' );

    for chaserIter = 1:MCsampleNum
        thisChaserLVLHPoint = QmatECItoLVLH*(thisPointVector(chaserIter,:)'-thisAvrgPoint');
        chaserTrajNormDevX(chaserIter,timeIter) = thisChaserLVLHPoint(1);
        chaserTrajNormDevY(chaserIter,timeIter) = thisChaserLVLHPoint(2);
        chaserTrajNormDevZ(chaserIter,timeIter) = thisChaserLVLHPoint(3);
    end
end

%%

figure(81)
set(gca,'FontSize',20)
set(gcf,'renderer','Painters')
hold on
grid on
title({'Experiment 2 Part 1:','Deviations in X-axis from Average Trajectory'})
for chaserIter = 1:MCsampleNum
    plot(1:length(MC_1_HPOP_ECI_X_Trajectories),chaserTrajNormDevX(chaserIter,:));
end
xlabel('Time [s]')
ylabel('Deviation [km]')
hold off

figure(82)
set(gca,'FontSize',20)
set(gcf,'renderer','Painters')
hold on
grid on
title({'Experiment 2 Part 1:','Deviations in Y-axis from Average Trajectory'})
for chaserIter = 1:MCsampleNum
    plot(1:length(MC_1_HPOP_ECI_Y_Trajectories),chaserTrajNormDevY(chaserIter,:));
end
xlabel('Time [s]')
ylabel('Deviation [km]')
hold off

figure(83)
set(gca,'FontSize',20)
set(gcf,'renderer','Painters')
hold on
grid on
title({'Experiment 2 Part 1:','Deviations in Z-axis from Average Trajectory'})
for chaserIter = 1:MCsampleNum
    plot(1:length(MC_1_HPOP_ECI_Z_Trajectories),chaserTrajNormDevZ(chaserIter,:));
end
xlabel('Time [s]')
ylabel('Deviation [km]')
hold off

%%

