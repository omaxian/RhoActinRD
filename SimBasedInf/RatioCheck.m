rng(0);
Noise = [0 1e-2];
PCASizes=40*ones(length(Noise),1);
RegZations=0*ones(length(Noise),1);

j=1;
ParInds = [3 6 10 11];
load('SyntheticXCors_ManySamps.mat')
load('PrincipalComponentsXCor.mat')
TestXCors=XCors_StarSim;
TestXCorsTr=U'*TestXCors;
TestXCorsTr=TestXCorsTr(1:max(PCASizes),:);
%bedge = -5:15;
ctrBin=round(mean(TestXCorsTr(j,:)));
bedge=ctrBin-4:0.5:ctrBin+4;
xpl = 1/2*(bedge(1:end-1)+bedge(2:end));
True=histcounts(TestXCorsTr(j,:),bedge);
True=True/((bedge(2)-bedge(1))*size(TestXCorsTr,2));
plot(xpl,True,'-k')
hold on
set(gca,'Colororderindex',1)

load('ScanUStimInPrior.mat')
nTrial=25;
nTest=1000;
testindsAll = randperm(length(AllParameters),nTest*nTrial);
AllCounts = zeros(length(Noise),length(xpl),nTrial);
for jEncSize=1:length(Noise)
for iTrial=1:nTrial
testinds = testindsAll((iTrial-1)*nTest+1:iTrial*nTest);
pChk = AllParameters(ParInds,testinds);
TestCors=U'*AllXCors(:,testinds);
TestCors=TestCors(1:40,:);
BinNum = 1+floor((TestCors(j,:)-bedge(1))/(bedge(2)-bedge(1)));
BinNum(BinNum<1)=1;
BinNum(BinNum>length(xpl))=length(xpl);
InputVec = [TestXCorsTr(:,1)'.*ones(nTest,1) pChk'];

load(strcat('TC_XCorNoise',num2str(Noise(jEncSize)),...
    'PCA',num2str(PCASizes(jEncSize)),...
    'Lam',num2str(RegZations(jEncSize)),'_Hyb.mat'),...
    'trainedClassifier')

[yfit,scores] = trainedClassifier.predictFcn(InputVec); 
LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
% Put the weights into the bins
TotalWeight = zeros(1,length(xpl));
for k=1:nTest    
    TotalWeight(BinNum(k))=TotalWeight(BinNum(k))+LikelihoodToEvidence(k);
end
AllCounts(jEncSize,:,iTrial)=TotalWeight;
AllCounts(jEncSize,:,iTrial)=AllCounts(jEncSize,:,iTrial)/((bedge(2)-bedge(1))...
    *sum(TotalWeight));

end
errorbar(xpl,mean(AllCounts(jEncSize,:,:),3),...
    2*std(AllCounts(jEncSize,:,:),[],3)/sqrt(nTrial),'LineWidth',2.0)
end



