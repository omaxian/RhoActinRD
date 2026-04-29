% This visualizes the weights and biases from a particular classifier.
Noise = [0 2e-3 5e-3 1e-2];
RegZations = [0 0 0 0];
PCASizes = 40*ones(1,length(Noise));
load('TC_XCorNoise0PCA40Lam0_Hyb5.mat')
WNone=trainedClassifier.ClassificationObj.LayerWeights{1,1};
[c,ed]=histcounts(WNone);
plot(1/2*(ed(1:end-1)+ed(2:end)),c)
hold on
for jEncSize=2:length(PCASizes)
    load(strcat('TC_XCorNoise',num2str(Noise(jEncSize)),'PCA',num2str(PCASizes(jEncSize)),...
        'Lam',num2str(RegZations(jEncSize)),'_Hyb5.mat'),'trainedClassifier')
    WOne{jEncSize-1}=trainedClassifier.ClassificationObj.LayerWeights{1,1};
    [c1]=histcounts(WOne{jEncSize-1},ed);
    plot(1/2*(ed(1:end-1)+ed(2:end)),c1)
end

load('TC_XCorNoise0.005PCA40Lam1e-05_Hyb5.mat')
WNone=trainedClassifier.ClassificationObj.LayerWeights{1,1};
c=histcounts(WNone,ed);
plot(1/2*(ed(1:end-1)+ed(2:end)),c,':k')

