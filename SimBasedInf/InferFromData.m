% This is the main file to infer the posterior from the data. 
% Assumes you already have a pre-trained classifier. 
% It will generate a plot of the marginal posterior over the parameters 
% you specify. Right now it is set up to do inference on the starfish data,
% holding kGAP constant
nP = 5;
figure(1);
PriorCutoff=0.5;
%tiledlayout(2,6,'Padding', 'none', 'TileSpacing', 'compact');
nModes = [10];
Noise =[0];
RegZations = [1e-3];
nLay = [20 40 20];
nW = 1;
for iDose=1
rng(1);
DataType='R';
ParInds=[2 3 6 10 11]; % (change to [2 3 6 10 11] to also vary kGAP)
FwdModel = @(p) RhoAndActinBasalNuc(p,1,0);
FwdModelPlot = @(p) RhoAndActinBasalNuc(p,1,1);
nPTot=11;
load('PrincipalComponentsXCor.mat')

if (DataType=='S')
    %ParamsTest = [0.7 0.4 30 1 1 0.817 0.055 0.65]; % wave
    %ParamsTest = [0.7 0.4 2 1 1 6 0.02 0.8]; % pulse
    %ParamsTest = [0.7 0.4 24.8508 1 1 0.8762 0.0536 0.8874]; % starfish ABC
    % ParamsTest = [0.7 0.4 2.14 1 1 4.69 0.0541 0.8519]; % worm ABC
    % 
    % figure(1);
    % tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
    % StatsTrue= FwdModelPlot(ParamsTest);
    % pbaspect([1 1 1])
    % title('True','interpreter','tex','FontWeight','Normal')
    % xticks(0:5:20)
    % xticklabels(xticks)
    % xlabel('$x$ ($\mu$m)')
    % ylabel('$t$ (s)')
    % TestXCors =StatsTrue.XCor(:);
    load('AllRepeats5.mat')
    ParamsTest=Allps(iDose,:);
    TestXCors = mean(XCorsPars(:,:,iDose),2);
elseif (DataType=='R')
    load('XCor_Starfish.mat')
    ParamsTest = [0.7 0.4 24.8508 1 1 0.8762 0.0536 0.8874]; % starfish ABC
    %load('XCor_Worm.mat')
    %ParamsTest = [0.7 0.4 2.14 1 1 4.69 0.0541 0.8519]; % worm ABC
    %load('XCorsByDose_Frog.mat')
    %DataXCor=XCorsByDose(:,:,iDose);
    TestXCors=DataXCor(:);
    %ParamsTest=[];
    % nexttile
    % imagesc(0:0.5:10,-120:5:120,DataXCor)
    % colormap(gca,'turbo')
    % clim([-1 1])
    % pbaspect([1 1 1])
    % xlabel('$\Delta r$ ($\mu$m)','interpreter','latex')
    % ylabel('$\Delta t$ (s)','interpreter','latex')
    % title('Data')
end
fixkgapval=ParamsTest(2);
try
ParamsTest=AddParams(ParamsTest);
catch
end

if (fixkgapval>0)
    EveryParameter = UniformParameterSets(ParInds(2:end),fixkgapval);
else
    EveryParameter = UniformParameterSets(ParInds);
end
% Pass through the prior
load(strcat('TCPriorOne_',num2str(nP),'.mat'))
[~,PriorScore]=TC_Sustainable.predictFcn(EveryParameter(:,ParInds));
InPrior = PriorScore(:,2) > PriorCutoff;
InPrior(EveryParameter(:,10)>1.6./EveryParameter(:,2))=0;
InPrior=logical(InPrior);
ParametersInPrior = EveryParameter(InPrior,:);

for jEncSize=1:length(nModes)
TestXCorsTr=U'*TestXCors;
TestXCorsTr = TestXCorsTr(1:nModes(jEncSize),:);
LikelihoodToEvidence = zeros(length(EveryParameter),1);
nSeed=2;
for seed=1:nSeed
load(strcat('NewTC',num2str(seed-1),'_XCorNoise',num2str(Noise),...
    'PCA',num2str(nModes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
% Use classifier to evaluate likelihood to evidence
InputVec = [TestXCorsTr'.*ones(size(ParametersInPrior,1),1) ParametersInPrior(:,ParInds)];
[~,scores] = trainedClassifier.predictFcn(InputVec); 
LikelihoodToEvidence(InPrior) = LikelihoodToEvidence(InPrior)+...
    1/nSeed*scores(:,2)./(1-scores(:,2));
if (DataType=='S')
   [yy,sc] = trainedClassifier.predictFcn(...
       [TestXCorsTr' ParamsTest(ParInds)])
end
end
[vals,inds]=sort(LikelihoodToEvidence);
mval = max(LikelihoodToEvidence);
[~,maxind]=max(LikelihoodToEvidence-mval*(PriorScore(:,2)<0.75));
% Compute marginals
figure(1);
[uvals{iDose},AllMarg{iDose}] = Compute2DMarginals(LikelihoodToEvidence,EveryParameter,...
   ParInds,1,ParamsTest,maxind,jEncSize,3);
end
XCorAvg=zeros(121,21);
nnz=0;
for j=1:10
    Stats1=RhoAndActinBasalNuc(EveryParameter(maxind,:),j,0,1);
    if (Stats1.EnoughExcitation)
        XCorAvg=XCorAvg+1/2*Stats1.XCor;
        nnz=nnz+1;
    end
    if (nnz==2)
        break;
    end
end
if (nnz~=2)
    keyboard
end
nexttile
imagesc(Stats1.rSim,Stats1.tSim,XCorAvg)
colormap(gca,'turbo')
xlabel('$\Delta r$ ($\mu$m)','interpreter','latex')
ylabel('$\Delta t$ (s)','interpreter','latex')
pbaspect([1 1 1])
clim([-1 1])
hold off
end