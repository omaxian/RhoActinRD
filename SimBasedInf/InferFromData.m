figure(2);
tiledlayout(1,5,'Padding', 'none', 'TileSpacing', 'compact');
for iDose=1:5
rng(0);
DataType='R';
ParInds=[3 6 10 11];
FixedParSet=0; % 0 for none, 1 to fix nucleation, 2 to fix assembly
FwdModel = @(p) RhoAndActinBasalNuc(p,1,0);
FwdModelPlot = @(p) RhoAndActinBasalNuc(p,1,1);
nPTot=11;
load('PrincipalComponentsXCor.mat')

if (DataType=='S')
    %ParamsTest = [0.7 0.4 30 1 1 0.817 0.055 0.65]; % wave
    ParamsTest = [0.7 0.4 2 1 1 6 0.02 0.8]; % pulse
    %ParamsTest = [0.7 0.4 24.8508 1 1 0.8762 0.0536 0.8874]; % starfish ABC
    %ParamsTest = [0.7 0.4 2.14 1 1 4.69 0.0541 0.8519]; % worm ABC
    if (FixedParSet==2)
        ParamsTest = [0.7 0.4 35 1 1 0.25 0.263589768883114 0.879455988011083];
        %ParamsTest = [0.7 0.4 35 1 1 0.25 0.214 1.514];
        %ParamsTest = [0.7 0.4 35 1 1 0.25 0.172 2.5874];
        %ParamsTest = [0.7 0.4 35 1 1 0.25 0.0563 9.0141];
    end

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
    %load('XCor_Starfish.mat')
    %ParamsTest = [0.7 0.4 24.8508 1 1 0.8762 0.0536 0.8874]; % starfish ABC
    %load('XCor_Worm.mat')
    %ParamsTest = [0.7 0.4 2.14 1 1 4.69 0.0541 0.8519]; % worm ABC
    load('XCorsByDose_Frog.mat')
    DataXCor=XCorsByDose(:,:,iDose);
    %ParamsTest = [0.7 0.4 35 1 1 0.25 0.249 0.898];
    %ParamsTest=[0.7 0.4 35 1 1 0.25 0.0563 9.0141];
    TestXCors=DataXCor(:);
    ParamsTest=[];
end
try
ParamsTest=AddParams(ParamsTest);
catch
end

% Compute probabilities
PCASizes = [40];%
Noise =[1e-2];
RegZations = [0];

%figure(5);
%tiledlayout(length(Noise),2,...
%  'Padding', 'none', 'TileSpacing', 'compact');

ParametersInPrior = ParameterSetsInPrior(ParInds,FixedParSet);
for jEncSize=1:length(PCASizes)
TestXCorsTr=U'*TestXCors;
TestXCorsTr = TestXCorsTr(1:PCASizes(jEncSize),:);
load(strcat('TC_XCorNoise',num2str(Noise(jEncSize)),'PCA',num2str(PCASizes(jEncSize)),...
    'Lam',num2str(RegZations(jEncSize)),'_Hyb.mat'),'trainedClassifier')

% Use classifier to evaluate likelihood to evidence
InputVec = [TestXCorsTr'.*ones(size(ParametersInPrior,1),1) ParametersInPrior(:,ParInds)];
[~,scores] = trainedClassifier.predictFcn(InputVec); 
LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
%if (DataType=='S')
%    [yy,sc] = trainedClassifier.predictFcn(...
%        [TestXCorsTr' ParamsTest(ParInds)])
%end

% Compute marginals
figure(2);
[uvals{iDose},AllMarg{iDose}] = Compute2DMarginals(LikelihoodToEvidence,ParametersInPrior,...
    ParInds,1,ParamsTest,FixedParSet);%,iDose,length(Doses));

% Samples from log likelihood
if (0)
SampVals = [1e3 1 1e-3];
PLInds = zeros(length(SampVals),1);
% figure(2*jEncSize);
% tiledlayout(2,length(SampVals),'Padding', 'none', 'TileSpacing', 'compact')
for jS=1:length(SampVals)
    [~,bind]=min(abs(vals-SampVals(jS)));
    if (bind==1)
        bind=10;
    end
    BestSustains=0;
    kmax=-10;
    k=0;
    nexttile(jS)
    while (~BestSustains && k >= kmax && bind+k>1)
        StatsBest = FwdModel(ParamsScan(inds(bind+k),:));
        BestSustains=StatsBest.EnoughExcitation;
        k=k-1;
        hold off
    end
    % colormap(gca,sky)
    % title(strcat('LER = $10^{',...
    %     num2str(round(log10(LikelihoodToEvidence(inds(bind+k+1)))*10)/10),'}$'))
    % xticklabels(xticks)
    % xlabel('$x$ ($\mu$m)','interpreter','latex')
    % if (jS==1)
    %     ylabel('$t$ (s)','interpreter','latex')
    % else
    %     yticklabels('')
    % end
    % pbaspect([1 1 1])
    % ylim([0 300])
    nexttile(1+jS+(jEncSize-1)*4)
    imagesc(0:0.5:10,-120:5:120,StatsBest.XCor)
    colormap(gca,turbo)
    pbaspect([1 1 1])
    if (jS==1)
        ylabel('$\Delta t$ (s)','interpreter','latex')
    else
        yticklabels('')
    end
    xlabel('$\Delta r$ ($\mu$m)','interpreter','latex')
    clim([-1 1])
    title(strcat('LER = $10^{',...
        num2str(round(log10(LikelihoodToEvidence(inds(bind+k+1)))*10)/10),'}$'))
    PLInds(jS)=inds(bind+k+1);
end
end

if (0)
figure
% Test set where we know the true result of fwd model
load('ScanUStimInPrior.mat')
nTest=1000;
testinds = randperm(length(AllParameters),nTest);
Ers = AllXCors(:,testinds)-TestXCors(:);
pChk = AllParameters(ParInds,testinds);
Ers = sqrt(sum(Ers.*Ers))./sqrt(sum(TestXCors(:).*TestXCors(:)));
% Sort the errors and take some evenly spaced ones
InputVec = [TestXCorsTr'.*ones(nTest,1) pChk'];
[yfit,scores] = trainedClassifier.predictFcn(InputVec); 
LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
%nexttile
scatter(Ers,log10(LikelihoodToEvidence),20,AllParameters(3,testinds),'filled');
tc = log10(LikelihoodToEvidence) > -10;
r = corrcoef(Ers(tc),log10(LikelihoodToEvidence(tc)));
legend(strcat('$r=$',num2str(round(r(1,2)*100)/100)))
if (jEncSize==length(PCASizes))
    xlabel('Norm of XCor error','interpreter','tex')
else
    xticklabels('')
end
ylabel('Log likelihood','interpreter','tex')
%ylim([-10 5])
pbaspect([1 1 1])
end
end
end