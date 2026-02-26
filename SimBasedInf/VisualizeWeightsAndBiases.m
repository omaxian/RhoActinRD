Noise = [0 0 0 0];
RegZations = [0 3.5e-4 1e-3 3.5e-3];
PCASizes = 40*ones(1,length(Noise));
load('TC_XCorNoise0PCA40Lam0_Hyb.mat')
WNone=trainedClassifier.ClassificationObj.LayerWeights{1,1};
[c,ed]=histcounts(WNone);
plot(1/2*(ed(1:end-1)+ed(2:end)),c)
hold on
for jEncSize=2:length(PCASizes)
    load(strcat('TC_XCorNoise',num2str(Noise(jEncSize)),'PCA',num2str(PCASizes(jEncSize)),...
        'Lam',num2str(RegZations(jEncSize)),'_Hyb.mat'),'trainedClassifier')
    WOne{jEncSize-1}=trainedClassifier.ClassificationObj.LayerWeights{1,1};
    [c1]=histcounts(WOne{jEncSize-1},ed);
    plot(1/2*(ed(1:end-1)+ed(2:end)),c1)
end

