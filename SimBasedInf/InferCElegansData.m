rng(2);
%% Load and proces experimental data
loadInData=1;
if (loadInData)
Name='nmy';
MovieNum=2;
pxlSize = 0.1;
FrRate = 1;
nModesFilt = 20;
MovAvgTime = 20;
L = 20;
Rho=double(load(strcat(Name,'Rho_',num2str(MovieNum),'.mat')).RhoData);
Actin=double(load(strcat(Name,'Actin_',num2str(MovieNum),'.mat')).ActinData);
Info=load(strcat(Name,'Info_',num2str(MovieNum),'.mat'));
FName=Info.Info.Name;
FrTime = Info.Info.TimeInt;
MaxT = ceil(200/FrTime);
% Adjust so that each frame has same mean
GlobalMeanRho=mean(Rho(:));
GlobalMeanActin=mean(Actin(:));
[ny,nx,nFr]=size(Rho);
for iT=1:nFr
    Rho(:,:,iT)=Rho(:,:,iT)-mean(mean(Rho(:,:,iT)))+GlobalMeanRho;
    Actin(:,:,iT)=Actin(:,:,iT)-mean(mean(Actin(:,:,iT)))+GlobalMeanActin;
end
% Moving average in time 
Rho = movmean(Rho,MovAvgTime,3);
Actin = movmean(Actin,MovAvgTime,3);
Rho = Rho(:,:,2*MovAvgTime:2*MovAvgTime+MaxT);
Actin = Actin(:,:,2*MovAvgTime:2*MovAvgTime+MaxT);
Nx=L/pxlSize;
x=(0:Nx-1)*pxlSize;
y=(0:Nx-1)*pxlSize;
% Filter data
[AllRho,RhoHatMean] = FilterData(Rho,nModesFilt);
[AllActin,ActHatMean] = FilterData(Actin,nModesFilt);
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
% Full movie 
% for g=1:size(AllRho,3)
% tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
% nexttile
% imagesc(x,y,AllRho(:,:,g))
% pbaspect([1 1 1])
% colormap sky
% clim([min(AllRho(:)) max(AllRho(:))])
% nexttile
% imagesc(x,y,AllActin(:,:,g))
% pbaspect([1 1 1])
% colormap sky
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
load('TrainedClassif_RelMeanFour1D10_3P_20K.mat')
% Random sampling
ParInd = [5 6 10];
nParCheck=2000;
LikelihoodToEvidence = zeros(nParCheck,nData);
ParamsScan = zeros(nParCheck,10);
PBounds = [0.01 0.5; 0 0.5; 0.5 10];
for jP=1:nParCheck
    p=PBounds(:,1)+(PBounds(:,2)-PBounds(:,1)).*rand(size(PBounds,1),1);
    Nreg1D = round(20/p(3));
    Areg = 400/Nreg1D^2;
    ParamsScan(jP,:)=[0.45 0.45 0.3*p(1) 4*p(1) 1*p(1) p(2) 0.05 1 0.1 Areg];
    % Use classifier to evaluate likelihood to evidence
    InputVec = [DataMeanRhoHat(gInds,:)'./DataMeanRhoHat(1,:)' ...
       DataMeanActHat(gInds,:)'./DataMeanActHat(1,:)' ...
       ParamsScan(jP,ParInd).*ones(nData,length(ParInd))];
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