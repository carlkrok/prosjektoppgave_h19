
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

coeffLarge = coeff*1000;

plotHandles = zeros(6,1);

figure(80)
hold on
axis equal
grid on
title('Experiment 2 Part 2: Chaser final position in target LVLH frame')
a = scatter3(x1(:,1), x1(:,2), x1(:,3), 'b');
a.MarkerEdgeAlpha = 0.1;
a.MarkerFaceAlpha = 0.1;
plotHandles(1) = plot3([mean1(1), mean1(1)+coeffLarge(1,1)], [mean1(2), mean1(2)+coeffLarge(2,1)], [mean1(3), mean1(3)+coeffLarge(3,1)], 'g');
plotHandles(2) = plot3([mean1(1), mean1(1)+coeffLarge(1,2)], [mean1(2), mean1(2)+coeffLarge(2,2)], [mean1(3), mean1(3)+coeffLarge(3,2)], 'g');
plotHandles(3) = plot3([mean1(1), mean1(1)+coeffLarge(1,3)], [mean1(2), mean1(2)+coeffLarge(2,3)], [mean1(3), mean1(3)+coeffLarge(3,3)], 'g');
%plotHandles(4) = plot_gaussian_ellipsoid(mean1, cov1, 3);
%b = plot_gaussian_ellipsoid(mean1(1:2), cov1(1:2,1:2), 3);
b = plot_gaussian_2d_xy(mean1, cov1, 3);
c = plot_gaussian_2d_xz(mean1, cov1, 3);
d = plot_gaussian_2d_yz(mean1, cov1, 3);
legend([plotHandles(1:3); b; c; d],'PC1','PC2','PC3','3\sigma_{xy}', '3\sigma_{xz}', '3\sigma_{yz}')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
%a.MarkerFaceAlpha = 0.2;
hold off









