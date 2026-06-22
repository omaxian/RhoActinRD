% This makes Fig. S5 - checking the ratio test. 
% It compares the distribution of cross correlations over the same
% parameter set p(x|theta) with random samples of cross correlation,
% weighted by the ratio the classifier gives you, p(x)r(x|theta)
nP=5;
nPv=5;
nLay=80;
nW=1;
load('AllRepeats5.mat')
offset=0;
nModesPlot=10;
nPerPlot = 5;
nTrial = 4;
Noise = 0;
PCASizes=2;
RegZations=0;
TrueMean=zeros(1,nModesPlot);
InfMean=zeros(nTrial,nModesPlot);
TrueStd=zeros(1,nModesPlot);
InfStd=zeros(nTrial,nModesPlot);
clear True AllWts
figure;
tiledlayout(1,5,'Padding', 'none', 'TileSpacing', 'compact');
rng(0);

for indrep=1
pTrue=Allps(indrep,:);
AllXCors=XCorsPars(:,:,indrep);
load('PrincipalComponentsXCor.mat')

pTrue = AddParams(pTrue); % starfish ABC
if (nP==2)
    ParInds = [3 6];
elseif (nP==5)
    ParInds = [2 3 6 10 11];
end
TruePar = pTrue(ParInds);
TestXCors=AllXCors;
TestXCorsTr=U'*TestXCors;

% Load data and classifier
for iTri=1:nTrial
if (iTri<nTrial/2)
try
load(strcat('TC0_XCorNoise',num2str(Noise),...
    'PCA',num2str(PCASizes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
catch
load(strcat('TC_XCorNoise',num2str(Noise),...
    'PCA',num2str(PCASizes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
end
else
try
load(strcat('TC1_XCorNoise',num2str(Noise),...
    'PCA',num2str(PCASizes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
catch
load(strcat('TC_XCorNoise',num2str(Noise),...
    'PCA',num2str(PCASizes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
end
end

load(strcat('UniformResults/Scan',num2str(nPv),'OneInPrior.mat'))
% Choose samples
nToCheck=5000;
inds=randperm(size(AllParameters,2),nToCheck);
AllParameters=AllParameters(:,inds);
AllXCors=AllXCors(:,inds);
nSamp = size(AllXCors,2);

TestCors=U'*AllXCors;
nModes = PCASizes;
% Two inference problems: infer parameters using the mean of p(x|theta),
% and infer data probabilities with fixed parameters
% Problem 1: parameters 
InputVec_Par = [repmat(mean(TestXCorsTr(1:nModes,:)'),nSamp,1) AllParameters(ParInds,:)'];
if (size(TestXCorsTr,2)==1)
    InputVec_Par = [repmat(TestXCorsTr(1:nModes)',nSamp,1) AllParameters(ParInds,:)'];
end
[yfit,scores] = trainedClassifier.predictFcn(InputVec_Par); 
LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
[~,plorder]=sort(LikelihoodToEvidence);
if (size(TestXCorsTr,2)==1)
    [~,ind]=max(LikelihoodToEvidence);
    pTrue=AllParameters(:,ind);
    TruePar=pTrue(ParInds,:)';
end
if (iTri==1)
if (nP==2)
    scatter(AllParameters(ParInds(1),:),AllParameters(ParInds(2),:),12,log10(LikelihoodToEvidence),'filled')
    hold on
    if (size(TestXCorsTr,2)>1)
    scatter(pTrue(ParInds(1)),pTrue(ParInds(2)),50,'ko','filled')
    end
else
    % nexttile(offset+1)
    % scatter(AllParameters(ParInds(1),plorder),AllParameters(ParInds(5),plorder),...
    %     12,log10(LikelihoodToEvidence(plorder)),'filled')
    % hold on
    % if (size(TestXCorsTr,2)>1)
    % scatter(pTrue(ParInds(1)),pTrue(ParInds(5)),50,'ko','filled')
    % end
    % g=clim;
    % clim([-3 ceil(max(g))])
    % t=turbo(100);
    % c=t(30:85,:);
    % colormap(gca,c)
    % title('Parameter space')
    % xlabel('$k_\textrm{GAP}$','interpreter','latex')
    % ylabel('$f_\rho$','interpreter','latex')
    % pbaspect([1 1 1])
    % nexttile(offset+2)
    % scatter(AllParameters(ParInds(4),plorder),AllParameters(ParInds(5),plorder),...
    %     12,log10(LikelihoodToEvidence(plorder)),'filled')
    % hold on
    % if (size(TestXCorsTr,2)>1)
    % scatter(pTrue(ParInds(4)),pTrue(ParInds(5)),50,'ko','filled')
    % end
    % g=clim;
    % clim([-3 ceil(max(g))])
    % t=turbo(100);
    % c=t(30:85,:);
    % colormap(gca,c)
    % title('Parameter space')
    % xlabel('$f_b$','interpreter','latex')
    % ylabel('$f_\rho$','interpreter','latex')
    %pbaspect([1 1 1])
    nexttile(offset+1)
    scatter(AllParameters(ParInds(2),plorder),AllParameters(ParInds(3),plorder),...
        12,log10(LikelihoodToEvidence(plorder)),'filled')
    hold on
    if (size(TestXCorsTr,2)>1)
    scatter(pTrue(ParInds(2)),pTrue(ParInds(3)),50,'ko','filled')
    end
end
g=clim;
clim([-3 ceil(max(g))])
t=turbo(100);
c=t(30:85,:);
colormap(gca,c)
title('Parameter space')
xlabel('$T_\textrm{fil}$','interpreter','latex')
ylabel('$\ell_\textrm{max}$','interpreter','latex')
pbaspect([1 1 1])
if (size(TestXCorsTr,2)==1)
    % Assuming expt data
    nexttile(offset+3)
    scatter(TestCors(1,:),TestCors(2,:),12,log10(LikelihoodToEvidence),'filled')
    hold on
    scatter(TestXCorsTr(1,:),TestXCorsTr(2,:),60,'ko','filled')
else
    InputVec_Data = [TestCors(1:nModes,:)' repmat(TruePar,nSamp,1)];
    [yfit,scores] = trainedClassifier.predictFcn(InputVec_Data); 
    LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
    if (sum(isinf(LikelihoodToEvidence))>0)
        LikelihoodToEvidence(isinf(LikelihoodToEvidence))=1e3;
        keyboard
    end
    [~,plorder]=sort(LikelihoodToEvidence);
    nexttile(offset+2)
    scatter(TestCors(1,plorder),TestCors(2,plorder),12,...
        log10(LikelihoodToEvidence(plorder)),'filled')
    hold on
    scatter(TestXCorsTr(1,:),TestXCorsTr(2,:),10,'ko','filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
end
g=clim;
clim([-3 ceil(max(g))])
t=turbo(100);
c=t(30:85,:);
colormap(gca,c)
title('PCA space')
xlabel('$\sigma_1$','interpreter','latex')
ylabel('$\sigma_2$','interpreter','latex')
pbaspect([1 1 1])
end

% Histograms
for j=1:nModesPlot
ctrBin=round(mean(TestXCorsTr(j,:)));
dx=0.5;
bedge=ctrBin-5:dx:ctrBin+5;
xpl = 1/2*(bedge(1:end-1)+bedge(2:end));
True(j,:) =histcounts(TestXCorsTr(j,:),bedge);
True(j,:) =True(j,:)/(dx*size(TestXCorsTr,2));
% Compute true means (first moment)
TrueMean(j) = sum(True(j,:).*xpl)*dx;
TrueStd(j) = sqrt(sum((xpl-TrueMean(j)).^2.*True(j,:))*dx);

AllCounts = zeros(length(xpl),1);
BinNum = 1+floor((TestCors(j,:)-bedge(1))/(bedge(2)-bedge(1)));
BinNum(BinNum<1)=1;
BinNum(BinNum>length(xpl))=length(xpl);
% Put the weights into the bins
TotalWeight = zeros(1,length(xpl));
for k=1:nSamp    
    TotalWeight(BinNum(k))=TotalWeight(BinNum(k))+LikelihoodToEvidence(k);
end
TotalWeight=TotalWeight/(dx*sum(TotalWeight));

InfMean(iTri,j) = sum(TotalWeight.*xpl)*dx;
InfStd(iTri,j) = sqrt(sum((xpl-InfMean(iTri,j)).^2.*TotalWeight)*dx);
AllWts(j,:,iTri)=TotalWeight;
end
end

plcut=0;
for j=1:nModesPlot
ctrBin=round(mean(TestXCorsTr(j,:)));
dx=0.5;
bedge=ctrBin-5:dx:ctrBin+5;
xpl = 1/2*(bedge(1:end-1)+bedge(2:end));
nexttile(offset+3+floor((j-1)/nPerPlot))
plot(xpl(True(j,:)>plcut),True(j,True(j,:)>plcut)+j,'-k')
hold on
set(gca,'Colororderindex',2)
MeanWt = mean(AllWts(j,:,:),3);
ErWt = 2*std(AllWts(j,:,:),[],3)/sqrt(nTrial);
ErWt=ErWt(MeanWt>plcut);
xpl=xpl(MeanWt>plcut);
MeanWt=MeanWt(MeanWt>plcut);
plot(xpl,MeanWt+j)
cplot=get(gca,'ColorOrder');
hold on
fill([xpl, fliplr(xpl)], [j+MeanWt-ErWt fliplr(j+MeanWt+ErWt)],cplot(2,:),...
    'FaceAlpha', 0.2, 'linestyle', 'none');
pbaspect([1 1 1])
legend('$p(x|\theta)$','$p(x)r(x,\theta)$','interpreter','latex')
xlabel('$\sigma_j$','interpreter','latex')
ylabel('$j$','interpreter','latex')
end
nexttile(offset+3+floor(nModesPlot/nPerPlot))
RelMeanEr = abs(InfMean-TrueMean)./InfStd;
plot(1:nModesPlot,mean(RelMeanEr),'-o')
hold on
fill([1:nModesPlot, fliplr(1:nModesPlot)], ...
    [mean(RelMeanEr)-2*std(RelMeanEr)/sqrt(nTrial) ...
    fliplr(mean(RelMeanEr)+2*std(RelMeanEr)/sqrt(nTrial))],cplot(1,:),...
    'FaceAlpha', 0.2, 'linestyle', 'none');
hold on
RelStdEr = InfStd./TrueStd;
set(gca,'Colororderindex',2)
plot(1:nModesPlot,mean(RelStdEr),'-o')
hold on
fill([1:nModesPlot, fliplr(1:nModesPlot)], ...
    [mean(RelStdEr)-2*std(RelStdEr)/sqrt(nTrial) ...
    fliplr(mean(RelStdEr)+2*std(RelStdEr)/sqrt(nTrial))],cplot(2,:),...
    'FaceAlpha', 0.2, 'linestyle', 'none');
grid on
pbaspect([1 1 1])

% Classifier ROC curve
% Cumulative weights
if (0)
WeightsBySamp = cumsum(LikelihoodToEvidence);
WeightsBySamp = WeightsBySamp/max(WeightsBySamp);
% Draw samples
nSamp = size(TestXCorsTr,2);
nTest = nSamp;
SampsFromRatio = zeros(size(TestCors,1),nSamp);
for g=1:nTest
    u = rand;
    index = find(WeightsBySamp<u,1,'last');
    if (isempty(index))
        index=0;
    end
    SampsFromRatio(:,g)=TestCors(:,index+1);
    indices(g)=index+1;
end


% Train classifier to distinguish two distributions
Weights=ones(nSamp+nTest,1);
TestData=[SampsFromRatio';TestXCorsTr'];
labels = [zeros(nTest,1); ones(nSamp,1)];
useNN=1;
CNetwork = fitcnet(TestData,labels,'LayerSizes', 25,'IterationLimit', 1000);

% 1. Get prediction scores (probabilities) from the trained network
% Scores will be an N-by-K matrix where K is the number of classes
[~,scores] =  predict(CNetwork, TestData);

% 3. Compute ROC metrics across all classes
rocObj = rocmetrics(labels, scores,[0 1]);

% 4. Plot the multiclass ROC curves
nexttile(offset+8)
plot(rocObj,'LineWidth',2,ClassNames=trainedClassifier.ClassificationObj.ClassNames(2));
hold on
end
end