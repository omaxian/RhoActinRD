% Checking the ratio test
% It compares the distribution of cross correlations over the same
% parameter set p(x|theta) with random samples of cross correlation,
% weighted by the ratio the classifier gives you, p(x)r(x|theta)
nP=5;
nPv=5;
nLay=[20 40 20];
nW=1;
offset=0;
nModesPlot=10;
nPerPlot = 5;
nTrial = 4;
Noise = 0;
nModes=10;
RegZations=1e-3;
nToCheck=20000; % Number of samples each time
TrueMean=zeros(1,nModesPlot);
InfMean=zeros(nTrial,nModesPlot);
TrueStd=zeros(1,nModesPlot);
InfStd=zeros(nTrial,nModesPlot);
clear True AllWts
if (offset==0)
figure;
tiledlayout(4,5,'Padding', 'none', 'TileSpacing', 'compact');
end
rng(0);
load('AllRepeats5.mat')
load('PrincipalComponentsXCor.mat')

for indrep=1:4
offset=5*(indrep-1);
% Define true parameters 
pTrue=AddParams(Allps(indrep,:));
ScanXCors=XCorsPars(:,:,indrep);
if (nP==2)
    ParInds = [3 6];
elseif (nP==5)
    ParInds = [2 3 6 10 11];
end
TruePar = pTrue(ParInds);
TestXCorsTr=U'*ScanXCors;

% Load data and classifier
for iTri=1:nTrial
if (iTri<=nTrial/2)
load(strcat('NewTC0_XCorNoise',num2str(Noise),...
    'PCA',num2str(nModes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
else
load(strcat('NewTC1_XCorNoise',num2str(Noise),...
    'PCA',num2str(nModes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
end
if (iTri==1)
    % Load parameter sets to check
    load(strcat('UniformResults/Scan',num2str(nPv),'NewOneInPrior.mat'))
end
% Choose samples
inds=randperm(size(AllParameters,2),nToCheck);
ScanParams=AllParameters(:,inds);
ScanXCors=AllXCors(:,inds);
TestCors=U'*ScanXCors;

% Two inference problems: infer parameters using the mean of p(x|theta),
% and infer data probabilities with fixed parameters
% Problem 1: parameters with data fixed
InputVec_Par = [repmat(mean(TestXCorsTr(1:nModes,:)'),nToCheck,1) ...
    ScanParams(ParInds,:)'];
[~,scores] = trainedClassifier.predictFcn(InputVec_Par); 
LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
[~,plorder]=sort(LikelihoodToEvidence);
% Plot the scatter plot for first trial only
if (iTri==1)
if (nP==2)
    nexttile(offset+1)
    scatter(ScanParams(ParInds(1),:),ScanParams(ParInds(2),:),12,log10(LikelihoodToEvidence),'filled')
    hold on
    scatter(pTrue(ParInds(1)),pTrue(ParInds(2)),50,'ko','filled')
else
    % nexttile(offset+1)
    % scatter(AllParameters(ParInds(1),plorder),AllParameters(ParInds(5),plorder),...
    %     12,log10(LikelihoodToEvidence(plorder)),'filled')
    % hold on
    % scatter(pTrue(ParInds(1)),pTrue(ParInds(5)),50,'ko','filled')
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
    % scatter(pTrue(ParInds(4)),pTrue(ParInds(5)),50,'ko','filled')
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
    scatter(ScanParams(ParInds(2),plorder),ScanParams(ParInds(3),plorder),...
        12,log10(LikelihoodToEvidence(plorder)),'filled')
    hold on
    scatter(pTrue(ParInds(2)),pTrue(ParInds(3)),50,'ko','filled')
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
colorbar
end

% Problem 2: data with parameters fixed
InputVec_Data = [TestCors(1:nModes,:)' repmat(TruePar,nToCheck,1)];
[yfit,scores] = trainedClassifier.predictFcn(InputVec_Data); 
LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
if (sum(isinf(LikelihoodToEvidence))>0)
    LikelihoodToEvidence(isinf(LikelihoodToEvidence))=1e3;
    keyboard
end
if (iTri==1)
    % Make the plot in PCA space
    [~,plorder]=sort(LikelihoodToEvidence);
    nexttile(offset+2)
    scatter(TestCors(1,plorder),TestCors(2,plorder),12,...
        log10(LikelihoodToEvidence(plorder)),'filled')
    hold on
    scatter(TestXCorsTr(1,:),TestXCorsTr(2,:),10,'ko','filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
    g=clim;
    clim([-3 ceil(max(g))])
    t=turbo(100);
    c=t(30:85,:);
    colormap(gca,c)
    title('PCA space')
    xlabel('$\gamma_1$','interpreter','latex')
    ylabel('$\gamma_2$','interpreter','latex')
    pbaspect([1 1 1])
    colorbar
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
for k=1:nToCheck    
    TotalWeight(BinNum(k))=TotalWeight(BinNum(k))+LikelihoodToEvidence(k);
end
TotalWeight=TotalWeight/(dx*sum(TotalWeight));

InfMean(iTri,j) = sum(TotalWeight.*xpl)*dx;
InfStd(iTri,j) = sqrt(sum((xpl-InfMean(iTri,j)).^2.*TotalWeight)*dx);
AllWts(j,:,iTri)=TotalWeight;
end
end

% Summary plot
plcut=-1e-16;
KLdiv=zeros(nTrial,nModesPlot);
for j=1:nModesPlot
ctrBin=round(mean(TestXCorsTr(j,:)));
bedge=ctrBin-5:dx:ctrBin+5;
xpl = 1/2*(bedge(1:end-1)+bedge(2:end));
nexttile(offset+3+floor((j-1)/nPerPlot))
plot(xpl(True(j,:)>plcut),True(j,True(j,:)>plcut)+j,'-k')
hold on
set(gca,'Colororderindex',2)
MeanWt = mean(AllWts(j,:,:),3);
% KL divergence of means
for iTri=1:nTrial
    KLdiv(iTri,j) = dx*sum(True(j,:).*log(True(j,:)./AllWts(j,:,iTri)),'omitnan');
end
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
if (j==nModesPlot)
legend('$p(x|\theta)$','$p(x)r(x,\theta)$','interpreter','latex','numColumns',1)
end
xlabel('$\gamma_j$','interpreter','latex')
ylabel('$j$','interpreter','latex')

end
nexttile(offset+3+floor(nModesPlot/nPerPlot))
RelMeanEr = abs(InfMean-TrueMean)./InfStd;
plot(1:nModesPlot,mean(RelMeanEr,1),'-o')
hold on
fill([1:nModesPlot, fliplr(1:nModesPlot)], ...
    [mean(RelMeanEr,1)-2*std(RelMeanEr,[],1)/sqrt(nTrial) ...
    fliplr(mean(RelMeanEr,1)+2*std(RelMeanEr,[],1)/sqrt(nTrial))],cplot(1,:),...
    'FaceAlpha', 0.2, 'linestyle', 'none');
hold on
RelStdEr = InfStd./TrueStd;
set(gca,'Colororderindex',2)
plot(1:nModesPlot,mean(RelStdEr,1),'-o')
hold on
fill([1:nModesPlot, fliplr(1:nModesPlot)], ...
    [mean(RelStdEr,1)-2*std(RelStdEr,[],1)/sqrt(nTrial) ...
    fliplr(mean(RelStdEr,1)+2*std(RelStdEr,[],1)/sqrt(nTrial))],cplot(2,:),...
    'FaceAlpha', 0.2, 'linestyle', 'none');
grid on
pbaspect([1 1 1])
xlabel('$j$','interpreter','latex')
set(gca,'Colororderindex',3)
KLdiv(isinf(KLdiv(:)))=3;
plot(1:nModesPlot,mean(KLdiv,1),'-o')
hold on
fill([1:nModesPlot, fliplr(1:nModesPlot)], ...
    [mean(KLdiv,1)-2*std(KLdiv,[],1)/sqrt(nTrial) ...
    fliplr(mean(KLdiv,1)+2*std(KLdiv,[],1)/sqrt(nTrial))],cplot(3,:),...
    'FaceAlpha', 0.2, 'linestyle', 'none');
legend('$|\hat \mu - \mu|/\hat \sigma$','','$\hat \sigma/\sigma$','','KL div','interpreter','latex','numColumns',1)

% Classifier ROC curve
% Cumulative weights
if (0)
WeightsBySamp = cumsum(LikelihoodToEvidence);
WeightsBySamp = WeightsBySamp/max(WeightsBySamp);
% Draw samples
nToCheck = size(TestXCorsTr,2);
nTest = nToCheck;
SampsFromRatio = zeros(size(TestCors,1),nToCheck);
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
Weights=ones(nToCheck+nTest,1);
TestData=[SampsFromRatio';TestXCorsTr'];
labels = [zeros(nTest,1); ones(nToCheck,1)];
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