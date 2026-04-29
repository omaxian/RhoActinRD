% This is the main file to infer parameters from data. 
% Assumes you already have a pre-trained classifier. 
% It will generate a plot of the average LER over the parameters 
% you specify. Right now it is set up to do inference on the starfish data,
% holding kGAP constant
figure(3);
%tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
for iDose=1
rng(0);
DataType='R';
ParInds=[3 6 10 11]; % (change to [2 3 6 10 11] to also vary kGAP)
FwdModel = @(p) RhoAndActinBasalNuc(p,1,0,1);
FwdModelPlot = @(p) RhoAndActinBasalNuc(p,1,1,1);
nPTot=11;
load('PrincipalComponentsXCor.mat')

if (DataType=='S')
    %ParamsTest = [0.7 0.4 30 1 1 0.817 0.055 0.65]; % wave
    %ParamsTest = [0.7 0.4 2 1 1 6 0.02 0.8]; % pulse
    %ParamsTest = [0.7 0.4 24.8508 1 1 0.8762 0.0536 0.8874]; % starfish ABC
    %ParamsTest = [0.7 0.4 2.14 1 1 4.69 0.0541 0.8519]; % worm ABC

    figure(1);
    tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
    StatsTrue= FwdModelPlot(ParamsTest);
    pbaspect([1 1 1])
    title('True','interpreter','tex','FontWeight','Normal')
    xticks(0:5:20)
    xticklabels(xticks)
    xlabel('$x$ ($\mu$m)')
    ylabel('$t$ (s)')
    TestXCors =StatsTrue.XCor(:);
elseif (DataType=='R')
    load('XCor_Starfish.mat')
    ParamsTest = [0.7 0.4 24.8508 1 1 0.8762 0.0536 0.8874]; % starfish ABC
    %load('XCor_Worm.mat')
    %ParamsTest = [0.7 0.4 2.14 1 1 4.69 0.0541 0.8519]; % worm ABC
    %load('XCorsByDose_Frog.mat')
    %DataXCor=XCorsByDose(:,:,iDose);
    TestXCors=DataXCor(:);
    %ParamsTest=[];
end
try
ParamsTest=AddParams(ParamsTest);
catch
end

% Compute probabilities
PCASizes = [40];%
Noise =[5e-3];
RegZations = [1e-5];


[EveryParameter,InPrior] = ParameterSetsInPrior(ParInds);
ParametersInPrior = EveryParameter(InPrior,:);
for jEncSize=1:length(PCASizes)
TestXCorsTr=U'*TestXCors;
TestXCorsTr = TestXCorsTr(1:PCASizes(jEncSize),:);
load(strcat('TC_XCorNoise',num2str(Noise(jEncSize)),'PCA',num2str(PCASizes(jEncSize)),...
    'Lam',num2str(RegZations(jEncSize)),'_Hyb5.mat'),'trainedClassifier')
% Use classifier to evaluate likelihood to evidence
InputVec = [TestXCorsTr'.*ones(size(ParametersInPrior,1),1) ParametersInPrior(:,ParInds)];
if (ParInds(1)>2)
    InputVec = [InputVec(:,1:PCASizes(jEncSize),:) 0.4*ones(size(InputVec,1),1) InputVec(:,PCASizes(jEncSize)+1:end)];
end
[~,scores] = trainedClassifier.predictFcn(InputVec); 
LikelihoodToEvidence = zeros(length(EveryParameter),1);
LikelihoodToEvidence(InPrior) = scores(:,2)./(1-scores(:,2));
if (DataType=='S')
   [yy,sc] = trainedClassifier.predictFcn(...
       [TestXCorsTr' 0.4 ParamsTest(ParInds)])
end
[vals,inds]=sort(LikelihoodToEvidence);
% Compute marginals
figure(3);
[uvals{iDose},AllMarg{iDose}] = Compute2DMarginals(LikelihoodToEvidence,EveryParameter,...
   ParInds,1,ParamsTest,jEncSize,3);
end
end