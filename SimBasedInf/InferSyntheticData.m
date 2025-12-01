Hybrid=1;
if (Hybrid)
    FwdModel = @(p) RhoAndActinBasalNuc(p,1,0);
    FwdModelPlot = @(p) RhoAndActinBasalNuc(p,1,1);
    ParamsTest = [0.7 0.4 26.9 1 1 0.817 0.055 0.65]; % starfish
    %ParamsTest = [0.7 0.4 4.9 1 1 1.6 0.202 1.072];
    ParamsTest=AddParams(ParamsTest);
    ParInds=[3 6 10 11];
    load('TC_XCor0_UStimInducHyb.mat')
    %load('TC_Four_wInducHyb.mat')
    load('TC_uStimInducSustain.mat')
    nPTot=11;
else
    FwdModel = @(p) RhoAndActinDiffusionPDEs(p,0.1,1,1,1,0);
    FwdModelPlot = @(p) RhoAndActinDiffusionPDEs(p,0.1,1,1,1,1); 
    ParamsTest = [0.7 0.4 26.9 1 1 0.817 0.055 0.65]; % starfish
    %ParamsTest = [0.7 0.4 4.9 1 1 1.6 0.202 1.072];
    ParInds = [1 5 6 10];
    load('TC_XCor_Cont.mat')
    %load('TC_Four_wInducHyb.mat')
    load('TC_ContSustain.mat')
    nPTot=10;
end
figure(1);
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
StatsTrue= FwdModelPlot(ParamsTest);
title('True')
TestMeanRhoHat = StatsTrue.MeanRhoHat(:);
TestMeanActHat = StatsTrue.MeanActinHat(:);
TestACorsRho = StatsTrue.ACorsRho(:);
TestACorsAct  = StatsTrue.ACorsAct(:);
TestXCors =StatsTrue.XCor(:);

% Using the trained classifier
nParCheck=5000;
ParamsScan=zeros(nParCheck,nPTot);
for iP=1:nParCheck
    InPrior=0;
    while (~InPrior)
        if (Hybrid)
            xr=rand(5,1);
            Params = [0.7; 0.4; 1+xr(2)*39; 1; 1; ...
                0.1+xr(3)*9.9; xr(4)*8; xr(5)*100]';
            MonomerClock = Params(6)/Params(4)+Params(6)/Params(5)+Params(3);
            NucRate_B = Params(7)/(MonomerClock*Params(6));
            NucRate_Rho = Params(8)/(MonomerClock*Params(6));
            Params(7)=NucRate_B;
            Params(8)=NucRate_Rho;
            Params=AddParams(Params);
        else
            PBounds = [0.45 0.75; 0.01 0.5; 0 0.5; 0.5 7];
            p=PBounds(:,1)+(PBounds(:,2)-PBounds(:,1)).*rand(size(PBounds,1),1);
            Nreg1D = round(20/p(4));
            Areg = 400/Nreg1D^2;
            Params=[p(1) 0.35 0.3*p(2) 4*p(2) p(2) p(3) 0.05 1 0.1 Areg];
        end
        InPrior = TC_Sustainable.predictFcn(Params(ParInds));
    end
    ParamsScan(iP,:)=Params;
end
% Use classifier to evaluate likelihood to evidence
Df = [TestMeanRhoHat(gInds)./TestMeanRhoHat(1); ...
    TestMeanActHat(gInds)./TestMeanActHat(1); ...
    TestACorsRho(cInds); ...
    TestACorsAct(cInds)];
InputVec = [Df'.*ones(nParCheck,1) ParamsScan(:,ParInds)];
%InputVec = [TestXCors(xInds,:)'.*ones(nParCheck,1) ParamsScan(:,ParInds)];
[yfit,scores] = trainedClassifier.predictFcn(InputVec); 
% [yy,sc] = trainedClassifier.predictFcn(...
%     [TestXCors(xInds,:)' ParamsTest(ParInds)])
LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
[vals,inds]=sort(LikelihoodToEvidence,'ascend');

% Statistics for ranked percentiles
BestSustains=0;
k=-1;
while (~BestSustains)
    k=k+1;
    BestInd=inds(end-k);
    StatsBest = FwdModel(ParamsScan(BestInd,:));
    BestSustains=StatsBest.EnoughExcitation;
end
% Visualize max likelihood trajectory
nexttile
StatsBest = FwdModelPlot(ParamsScan(BestInd,:));
title('Most likely')

% 90th prctile
NinetySustains=0;
k=-1;
while (~NinetySustains)
    k=k+1;
    NinetyInd=inds(end-nParCheck/10-k);
    StatsNinety = FwdModel(ParamsScan(NinetyInd,:));
    NinetySustains=StatsNinety.EnoughExcitation;
end

% 75th prctile
SFSustains=0;
k=-1;
while (~SFSustains)
    k=k+1;
    SevFInd=inds(end-nParCheck/4-k);
    StatsSF = FwdModel(ParamsScan(SevFInd,:));
    SFSustains=StatsSF.EnoughExcitation;
end



%Plot clusters in parameter space 
figure(2);
tiledlayout(4,2,'Padding', 'none', 'TileSpacing', 'compact');
if (Hybrid)
Tilex = [6 10];
Tiley = [3 11];
xLabels = ["$\ell$" "$f_b$"];
yLabels = ["$T_\textrm{fil}$" "$f_\rho$"];
else
Tilex = [1 6];
Tiley = [5 10];
xLabels = ["$k_\textrm{off}^{(0)}$" "$D_f$"];
yLabels = ["$k_\textrm{diss}$" "$A_\textrm{reg}$"];
end
Cutoff = max(log10(LikelihoodToEvidence))-10;
Goodinds = inds(log10(LikelihoodToEvidence(inds))>Cutoff);
for iTile=1:2
nexttile
scatter(ParamsScan(Goodinds,Tilex(iTile)),ParamsScan(Goodinds,Tiley(iTile)),10,...
    log10(LikelihoodToEvidence(Goodinds)),'filled')
hold on
scatter(ParamsTest(Tilex(iTile)),ParamsTest(Tiley(iTile)),...
    100,'ks','filled')
scatter(ParamsScan(BestInd,Tilex(iTile)),ParamsScan(BestInd,Tiley(iTile)),...
    100,'m^','filled')
ylabel(yLabels(iTile))
xlabel(xLabels(iTile))
colormap jet
clim([Cutoff max(log10(LikelihoodToEvidence))])
pbaspect([1 1 1])
end
colorbar

% Statistics
nM=10;
gInds=(1:nM)+(0:10:(nM-1)*10)';
cInds=gInds(1:10)+200;
gInds=gInds(2:10);
nexttile
plot(1:nM-1,StatsTrue.MeanRhoHat(gInds(:))./StatsTrue.MeanRhoHat(1),'-k')
hold on
plot(1:nM-1,StatsBest.MeanRhoHat(gInds(:))./StatsBest.MeanRhoHat(1),'Color',[0 0.44 0.89 1])
plot(1:nM-1,StatsNinety.MeanRhoHat(gInds(:))./StatsNinety.MeanRhoHat(1),'Color',[0 0.44 0.89 0.5])
plot(1:nM-1,StatsSF.MeanRhoHat(gInds(:))./StatsSF.MeanRhoHat(1),'Color',[0 0.44 0.89 0.25])
pbaspect([1 1 1])
xlabel('$m$')
ylabel('$<|\hat \rho(m,0)|>$')
nexttile
plot(1:nM-1,StatsTrue.MeanActinHat(gInds(:))./StatsTrue.MeanActinHat(1),'-k')
hold on
plot(1:nM-1,StatsBest.MeanActinHat(gInds(:))./StatsBest.MeanActinHat(1),'Color',[0 0.44 0.89 1])
plot(1:nM-1,StatsNinety.MeanActinHat(gInds(:))./StatsNinety.MeanActinHat(1),'Color',[0 0.44 0.89 0.5])
plot(1:nM-1,StatsSF.MeanActinHat(gInds(:))./StatsSF.MeanActinHat(1),'Color',[0 0.44 0.89 0.25])
xlabel('$m$')
ylabel('$<|\hat f(m,0)|>$')
pbaspect([1 1 1])
nexttile
plot(0:nM-1,StatsTrue.ACorsRho(cInds(:)),'-k')
hold on
plot(0:nM-1,StatsBest.ACorsRho(cInds(:)),'Color',[0 0.44 0.89 1])
plot(0:nM-1,StatsNinety.ACorsRho(cInds(:)),'Color',[0 0.44 0.89 0.5])
plot(0:nM-1,StatsSF.ACorsRho(cInds(:)),'Color',[0 0.44 0.89 0.25])
xlabel('$m$')
ylabel('acor($|\hat \rho|,\Delta t = 4$)')
pbaspect([1 1 1])
nexttile
plot(0:nM-1,StatsTrue.ACorsAct(cInds(:)),'-k')
hold on
plot(0:nM-1,StatsBest.ACorsAct(cInds(:)),'Color',[0 0.44 0.89 1])
plot(0:nM-1,StatsNinety.ACorsAct(cInds(:)),'Color',[0 0.44 0.89 0.5])
plot(0:nM-1,StatsSF.ACorsAct(cInds(:)),'Color',[0 0.44 0.89 0.25])
xlabel('$m$')
ylabel('acor($|\hat f|,\Delta t = 4$)')
pbaspect([1 1 1])
nexttile
ResampledT = -120:2:120;
ResampledX = 0:0.5:10;
[X,T]=meshgrid(ResampledX,ResampledT);
X0Inds=find(X==0);
T0Inds=find(T==0);
plot(T(X0Inds),StatsTrue.XCor(X0Inds),'-k')
hold on
plot(T(X0Inds),StatsBest.XCor(X0Inds),'Color',[0 0.44 0.89 1])
plot(T(X0Inds),StatsNinety.XCor(X0Inds),'Color',[0 0.44 0.89 0.5])
plot(T(X0Inds),StatsSF.XCor(X0Inds),'Color',[0 0.44 0.89 0.25])
xlabel('$\Delta t$')
ylabel('$R_{\rho f}(0,\Delta t)$')
pbaspect([1 1 1])
nexttile
plot(X(T0Inds),StatsTrue.XCor(T0Inds),'-k')
hold on
plot(X(T0Inds),StatsBest.XCor(T0Inds),'Color',[0 0.44 0.89 1])
plot(X(T0Inds),StatsNinety.XCor(T0Inds),'Color',[0 0.44 0.89 0.5])
plot(X(T0Inds),StatsSF.XCor(T0Inds),'Color',[0 0.44 0.89 0.25])
xlabel('$\Delta r$')
ylabel('$R_{\rho f}(\Delta r,0)$')
pbaspect([1 1 1])