% % Process the data into a list
% % Make a list of 2D cross correlations and parameters
% % PRE-ALLOCATE!
% nSampTot=20000;
% AllXCors=zeros(2541,nSampTot);
% AllMeanRhoHat=zeros(100,nSampTot);
% AllMeanActHat=zeros(100,nSampTot);
% AllACorsRho=zeros(400,nSampTot);
% AllACorsAct=zeros(400,nSampTot);
% AllParams=zeros(nSampTot,10);
% %PctAboveSaddle = zeros(2001,nSampTot);
% for jS=1:100
%     %try
%     load(strcat('Scan3PDbl_',num2str(jS),'.mat'));
%     for k=1:nSamp
%         AllXCors(:,(jS-1)*nSamp+k) = AllStats{k}.XCor(:);
%         AllMeanRhoHat(:,(jS-1)*nSamp+k) = AllStats{k}.MeanRhoHat(:);
%         AllMeanActHat(:,(jS-1)*nSamp+k) = AllStats{k}.MeanActinHat(:);
%         AllACorsRho(:,(jS-1)*nSamp+k) = AllStats{k}.ACorsRho(:);
%         AllACorsAct(:,(jS-1)*nSamp+k) = AllStats{k}.ACorsAct(:);
%         PctAboveSaddle{(jS-1)*nSamp+k} = AllStats{k}.PctAboveSaddle;
%         AllParams((jS-1)*nSamp+k,:) = Params(k,:);
%     end
%    % catch
%    %     jS
%    % end
% end
% %[~,inds]=unique(AllParams,'rows');
% % AllParams=AllParams(inds,:);
% % PctAboveSaddle=PctAboveSaddle(:,inds);
% % AllMeanRhoHat=AllMeanRhoHat(:,inds);
% % AllMeanActHat=AllMeanActHat(:,inds);
% % AllACorsRho=AllACorsRho(:,inds);
% % AllACorsAct=AllACorsAct(:,inds);
% % AllXCors=AllXCors(:,inds);
% clearvars -except PctAboveSaddle AllParams AllMeanRhoHat AllMeanActHat AllACorsRho AllACorsAct AllXCors
% save('DataForClassif_ThreeParamsk45Dbl.mat')
% return
% 
% %% Generate fake data and compare against parameter sets
% % How many of the parameter sets are just sitting in the steady state?
% % Remove parameter sets that are doing nothing
% % % 
%Try to identify the minima again
%Generate some data for testing
if (1)
Dfs=[0.01 0.1 0.4];
koffs=[0.02 0.1 0.4];
Lreg = [10/3];
[Dfs,koffs,Lreg]=meshgrid(Dfs,koffs,Lreg);
Dfs=Dfs(:);
koffs=koffs(:);
Lreg=Lreg(:);
%Lreg=[10/3 1 5 2 10/3 10 2.5 2.5 10/6]';
%Dfs=Dfs([1 3 6 7 9]);
%koffs=koffs([1 3 6 7 9]);
%Lreg=Lreg([1 3 6 7 9]);
nData = length(Lreg);
nP=10;
ParamsTest = zeros(nData,nP);
figure;
tiledlayout(3,3,...
    'Padding', 'none', 'TileSpacing', 'compact')
clear TestMeanRhoHat TestMeanActHat TestACorsRho TestACorsAct TestXCors
for iD=1:nData
    p=[0.45 0.45 0.3 4 1 0.5 0.05 1 0.1 10];
    p(3:5)=[0.3 4 1]*koffs(iD);
    p(6)=Dfs(iD);
    p(10)=Lreg(iD)^2;
    ParamsTest(iD,:)=p;
    Stats=RhoAndActinPDEMod(ParamsTest(iD,:),0.1,1,1,2,1);
    TestMeanRhoHat(:,iD) = Stats.MeanRhoHat(:);
    TestMeanActHat(:,iD) = Stats.MeanActinHat(:);
    TestACorsRho(:,iD) = Stats.ACorsRho(:);
    TestACorsAct(:,iD) = Stats.ACorsAct(:);
    TestXCors(:,iD)=Stats.XCor(:);
    if (mean(Stats.PctAboveSaddle)>0.99 || mean(Stats.PctAboveSaddle)<0.01)
        TestXCors(:,iD)=TestXCors(:,iD)*0;
        TestACorsRho(:,iD)=TestACorsRho(:,iD)*0;
        TestACorsAct(:,iD)=TestACorsAct(:,iD)*0;
    end
end
end
% 
if (0)
% Plot the true difference in the Xcors
nSampOG=size(AllParams,1);
LikelihoodToEvidence = zeros(nSampOG,nData);
for jChk=1:nData
    for iSamp=1:nSampOG
        LikelihoodToEvidence(iSamp,jChk)=...
            10.^(-norm(AllXCors(:,iSamp)-TestXCors(:,jChk)));
    end
end
nToChk = 5000;
p = randperm(nSampOG,nToChk);
figure
ScatterPlotParams(LikelihoodToEvidence(p,:),AllParams(p,:),...
    ParamsTest);
end


% Using the trained classifier
load('TrainedClassif_RelMeanFour1D10_2P.mat')
% Random sampling
ParInd = [5 6];
nParCheck=2000;
LikelihoodToEvidence = zeros(nParCheck,nData);
ParamsScan = zeros(nParCheck,nP);
PBounds = [0.01 0.5; 0 0.5; 0.5 10];
for jP=1:nParCheck
    p=PBounds(:,1)+(PBounds(:,2)-PBounds(:,1)).*rand(size(PBounds,1),1);
    Nreg1D = round(20/p(3));
    Areg = 400/Nreg1D^2;
    if (sum(ParInd==10)==0)
        Areg = (10/3)^2;
    end
    ParamsScan(jP,:)=[0.45 0.45 0.3*p(1) 4*p(1) 1*p(1) p(2) 0.05 1 0.1 Areg];
    % Use classifier to evaluate likelihood to evidence
    %InputVec = [TestXCors(xInds,:)' ParamsScan(jP,ParInd).*ones(nData,length(ParInd))];
    InputVec = [TestMeanRhoHat(gInds,:)'./TestMeanRhoHat(1,:)' ...
       TestMeanActHat(gInds,:)'./TestMeanActHat(1,:)' ...
       ParamsScan(jP,ParInd).*ones(nData,length(ParInd))];
    %InputVec = [TestACorsRho(cInds,:)' TestACorsAct(cInds,:)' ...
    % ParamsScan(jP,ParInd).*ones(nData,length(ParInd))];
    %InputVec = [TestMeanActHat(gInds,:)'./TestMeanActHat(1,:)' ...
    %   TestACorsAct(cInds,:)' ...
    %   ParamsScan(jP,ParInd).*ones(nData,length(ParInd))];
    [yfit,scores] = trainedClassifier.predictFcn(InputVec); 
    LikelihoodToEvidence(jP,:) = scores(:,2)./(1-scores(:,2));
end

%Identify most likely parameter set
%Plot clusters in parameter space 
figure
IndsChosen=ScatterPlotParams(LikelihoodToEvidence,ParamsScan,ParamsTest);
figure
tiledlayout(3,3,...
    'Padding', 'none', 'TileSpacing', 'compact')
for jP=1:nData
    ind=IndsChosen(jP);
    Stats = RhoAndActinPDEMod(ParamsScan(ind,:),0.1,1,1,1,1);
end
return
figure
tiledlayout(nData,2,...
    'Padding', 'none', 'TileSpacing', 'compact');
for iD=1:nData
[vals,inds]=sort(LikelihoodToEvidence(:,iD));
Stats22 = RhoAndActinPDEMod(ParamsTest(iD,:),0.1,1,1,2,0);
Stats0 = RhoAndActinPDEMod(ParamsScan(inds(end),:),0.1,1,1,1,0);
Stats1 = RhoAndActinPDEMod(ParamsScan(inds(end-100),:),0.1,1,1,1,0);
Stats2 = RhoAndActinPDEMod(ParamsScan(inds(end-250),:),0.1,1,1,1,0);
nexttile
plot(Stats22.MeanRhoHat(gInds)./Stats22.MeanRhoHat(1),'-k')
hold on
plot(Stats0.MeanRhoHat(gInds)./Stats0.MeanRhoHat(1),'Color',[0 0.44 0.89 1])
plot(Stats1.MeanRhoHat(gInds)./Stats1.MeanRhoHat(1),'Color',[0 0.44 0.89 0.5])
plot(Stats2.MeanRhoHat(gInds)./Stats2.MeanRhoHat(1),'Color',[0 0.44 0.89 0.25])
nexttile
plot(Stats22.MeanActinHat(gInds)./Stats22.MeanActinHat(1),'-k')
hold on
plot(Stats0.MeanActinHat(gInds)./Stats0.MeanActinHat(1),'Color',[0 0.44 0.89 1])
plot(Stats1.MeanActinHat(gInds)./Stats1.MeanActinHat(1),'Color',[0 0.44 0.89 0.5])
plot(Stats2.MeanActinHat(gInds)./Stats2.MeanActinHat(1),'Color',[0 0.44 0.89 0.25])
end