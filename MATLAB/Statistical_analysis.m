x1 = relEndPosLVLHHPOP_1_chaser.*10^3;
%x2 = absDeviationEndPosHPOP_1;


%%


mean1 = mean(x1)
%mean2 = mean(x2)

%%


var1 = var(x1)
%var2 = var(x2)

%%

skew1 = skewness(x1)
%skew2 = skewness(x2)

%%

kurt1 = kurtosis(x1)
%kurt2 = kurtosis(x2)

%% 

cov(relEndPosLVLHHPOP_1_chaser.*10^3)
