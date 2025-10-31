rng(0);
%% Load and proces experimental data
loadInData=1;
if (loadInData)
pxlSize = 5/9;
FrTime = 4;
nModesFilt = 20;
L = 20;
Rho=tiffreadVolume("BementRho.tif");
Actin=tiffreadVolume("BementRGA.tif");
Nx=L/pxlSize;
Rho=Rho(11:10+Nx,11:10+Nx,:);
Actin=Actin(11:10+Nx,11:10+Nx,:);
% Difference subtraction (removes static signal)
Rho=Rho-Rho(:,:,end);
Actin=Actin-Actin(:,:,end);
Rho=Rho(:,:,1:end-1);
Actin=Actin(:,:,1:end-1);
% Moving average in time 
% Rho = movmean(Rho,MovAvgTime,3);
% Actin = movmean(Actin,MovAvgTime,3);
% Rho = Rho(:,:,2*MovAvgTime:2*MovAvgTime+MaxT);
% Actin = Actin(:,:,2*MovAvgTime:2*MovAvgTime+MaxT);
x=(0:Nx-1)*pxlSize;
y=(0:Nx-1)*pxlSize;
% Filter data
[AllRho,RhoHatMean] = FilterData(Rho,nModesFilt);
[AllActin,ActHatMean] = FilterData(Actin,nModesFilt);
[UvalsF,dtvalsF,DistsByRF] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    AllRho,AllActin,1);
DistsByRF=DistsByRF/max(abs(DistsByRF(:)));
ResampledT = -120:2:120;
ResampledX = 0:0.5:10;
[X,T]=meshgrid(ResampledX,ResampledT);
X0Inds=find(X==0 & T>=-60 & T<=60);
T0Inds=find(T==0);
DataXCor=ResampleXCor(DistsByRF,dtvalsF,UvalsF,ResampledX,ResampledT,11,121);
% Make a kymograph
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
ax1=nexttile;
C2=[0.87 0.49 0];
C1=[0.95 0.9 0.9];
Cmap=C1+(0:100)'/100.*(C2-C1);
RhoKymo=reshape(AllRho(3*Nx/4,:,:),Nx,size(AllRho,3))';
ActKymo=reshape(AllActin(3*Nx/4,:,:),Nx,size(AllActin,3))';
imagesc((0:Nx-1)*pxlSize,(0:size(RhoKymo,1)-1)*FrTime,RhoKymo)
clim([min(AllRho(:)) max(AllRho(:))])
xlabel('$x$ ($\mu$m)')
ylabel('$t$ (s)')
pbaspect([1 1 1])
colormap(ax1,sky)
xticklabels('')
ax2=nexttile;
imagesc((0:Nx-1)*pxlSize,(0:size(RhoKymo,1)-1)*FrTime,ActKymo)
clim([min(AllActin(:)) max(AllActin(:))])
xlabel('$x$ ($\mu$m)')
%ylabel('$t$ (s)')
pbaspect([1 1 1])
yticklabels('')
xticklabels('')
colormap(ax2,Cmap)
%Full movie 
% for g=1:size(AllRho,3)
% tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
% ax3=nexttile;
% imagesc(x,y,AllRho(:,:,g))
% pbaspect([1 1 1])
% colormap(ax3, sky)
% clim([min(AllRho(:)) max(AllRho(:))])
% ax4=nexttile;
% imagesc(x,y,AllActin(:,:,g))
% pbaspect([1 1 1])
% colormap(ax4,Cmap)
% clim([min(AllActin(:)) max(AllActin(:))])
% drawnow
% end
nModesForClassif=10;
DataMeanRhoHat=RhoHatMean(1:nModesForClassif,1:nModesForClassif);
DataMeanRhoHat=DataMeanRhoHat(:);
DataMeanActHat=ActHatMean(1:nModesForClassif,1:nModesForClassif);
DataMeanActHat=DataMeanActHat(:);
end

% load('DataForClassif_ThreeParamsk45.mat')
% AllMeanRhos=AllMeanRhoHat(gInds,:)./AllMeanRhoHat(1,:);
% AllMeanActs=AllMeanActHat(gInds,:)./AllMeanActHat(1,:);
% DRho = DataMeanRhoHat(gInds)/DataMeanRhoHat(1);
% DAct = DataMeanActHat(gInds)/DataMeanActHat(1);
% diffnorm = [AllMeanRhos;AllMeanActs]-[DRho;DAct];

nData=1;
% Using the trained classifier
load('TrainedClassif_RelMeanFour1D10_2P.mat')
% Random sampling
ParInd = [5 6];
nParCheck=2000;
LikelihoodToEvidence = zeros(nParCheck,nData);
ParamsScan = zeros(nParCheck,10);
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
    InputVec = [DataMeanRhoHat(gInds,:)'./DataMeanRhoHat(1,:)' ...
       DataMeanActHat(gInds,:)'./DataMeanActHat(1,:)' ...
       ParamsScan(jP,ParInd).*ones(nData,length(ParInd))];
    %InputVec = [XCorData(xInds)' ParamsScan(jP,ParInd).*ones(nData,length(ParInd))];
    [yfit,scores] = trainedClassifier.predictFcn(InputVec); 
    LikelihoodToEvidence(jP,:) = scores(:,2)./(1-scores(:,2));
end

%Identify most likely parameter set
%Plot clusters in parameter space 
figure
IndsChosen=ScatterPlotParams(LikelihoodToEvidence,ParamsScan);
figure
tiledlayout(nData,2,...
    'Padding', 'none', 'TileSpacing', 'compact')
for jP=1:nData
    ind=IndsChosen(jP);
    Stats = RhoAndActinPDEMod(ParamsScan(ind,:),0.1,1,1,1,1);
end

figure
tiledlayout(nData,2,...
    'Padding', 'none', 'TileSpacing', 'compact');
for iD=1:nData
[vals,inds]=sort(LikelihoodToEvidence(:,iD));
Stats0 = RhoAndActinPDEMod(ParamsScan(inds(end),:),0.1,1,1,1,0);
Stats1 = RhoAndActinPDEMod(ParamsScan(inds(end-100),:),0.1,1,1,1,0);
Stats2 = RhoAndActinPDEMod(ParamsScan(inds(end-250),:),0.1,1,1,1,0);
% nexttile
% plot(T(X0Inds),DataXCor(X0Inds),'-k')
% hold on
% plot(T(X0Inds),Stats0.XCor(X0Inds),'Color',[0 0.44 0.89 1])
% plot(T(X0Inds),Stats1.XCor(X0Inds),'Color',[0 0.44 0.89 0.5])
% plot(T(X0Inds),Stats2.XCor(X0Inds),'Color',[0 0.44 0.89 0.25])
% xlabel('$\Delta t$')
% ylabel('$R_{\rho f}(0,\Delta t)$')
% nexttile
% plot(X(T0Inds),DataXCor(T0Inds),'-k')
% hold on
% plot(X(T0Inds),Stats0.XCor(T0Inds),'Color',[0 0.44 0.89 1])
% plot(X(T0Inds),Stats1.XCor(T0Inds),'Color',[0 0.44 0.89 0.5])
% plot(X(T0Inds),Stats2.XCor(T0Inds),'Color',[0 0.44 0.89 0.25])
% xlabel('$\Delta r$')
% ylabel('$R_{\rho f}(\Delta r,0)$')
nexttile
plot(DataMeanRhoHat(gInds)./DataMeanRhoHat(1),'-k')
hold on
plot(Stats0.MeanRhoHat(gInds)./Stats0.MeanRhoHat(1),'Color',[0 0.44 0.89 1])
plot(Stats1.MeanRhoHat(gInds)./Stats1.MeanRhoHat(1),'Color',[0 0.44 0.89 0.5])
plot(Stats2.MeanRhoHat(gInds)./Stats2.MeanRhoHat(1),'Color',[0 0.44 0.89 0.25])
nexttile
plot(DataMeanActHat(gInds)./DataMeanActHat(1),'-k')
hold on
plot(Stats0.MeanActinHat(gInds)./Stats0.MeanActinHat(1),'Color',[0 0.44 0.89 1])
plot(Stats1.MeanActinHat(gInds)./Stats1.MeanActinHat(1),'Color',[0 0.44 0.89 0.5])
plot(Stats2.MeanActinHat(gInds)./Stats2.MeanActinHat(1),'Color',[0 0.44 0.89 0.25])
end