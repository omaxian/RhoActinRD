% Process and filter experimental data for starfish and worms
rng(0);
%% Starfish
nModesFilt = 50;
pxlSize = 5/9;
FrTime = 4;
Rho=tiffreadVolume("BementRho.tif");
Actin=tiffreadVolume("BementRGA.tif");
Nx=size(Rho,1);
%Rho=Rho(11:10+Nx,11:10+Nx,:);
%Actin=Actin(11:10+Nx,11:10+Nx,:);
% % Difference subtraction (removes static signal)
Rho=Rho-Rho(:,:,end);
Actin=Actin-Actin(:,:,end);
Rho=Rho(:,:,1:end-1);
Actin=Actin(:,:,1:end-1);
% Adjust so that each frame has same mean
GlobalMeanRho=mean(Rho(:));
GlobalMeanActin=mean(Actin(:));
[ny,nx,nFr]=size(Rho);
for iT=1:nFr
    Rho(:,:,iT)=Rho(:,:,iT)-mean(mean(Rho(:,:,iT)))+GlobalMeanRho;
    Actin(:,:,iT)=Actin(:,:,iT)-mean(mean(Actin(:,:,iT)))+GlobalMeanActin;
end
Rho=double(Rho);
Actin=double(Actin);

%% Worms 
if (1)
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
MaxT = ceil(300/FrTime);
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
end

%% Filtrring and statistics
x=(0:Nx-1)*pxlSize;
y=(0:Nx-1)*pxlSize;
% Filter data
[AllRho,RhoHatMean,ACorsRho,TimeLags] = FilterData(Rho,nModesFilt,FrTime);
[AllActin,ActHatMean,ACorsAct,TimeLags] = FilterData(Actin,nModesFilt,FrTime);
[UvalsF,dtvalsF,DistsByRF] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    AllRho,AllActin,1);
DistsByRF=DistsByRF/max(abs(DistsByRF(:)));
ResampledT = -120:2:120;
ResampledX = 0:0.5:10;
[X,T]=meshgrid(ResampledX,ResampledT);
X0Inds=find(X==0);
T0Inds=find(T==0);
DataXCor=ResampleXCor(DistsByRF,dtvalsF,UvalsF,ResampledX,ResampledT,11,121);
% Make a kymograph
figure;
tiledlayout(3,1,'Padding', 'none', 'TileSpacing', 'compact');
ax0=nexttile;
imagesc((0:Nx-1)*pxlSize,(0:Nx-1)*pxlSize,AllRho(:,:,1))
colormap(ax0,sky)
xlabel('$x$ ($\mu$m)')
title('Snapshot','Interpreter','tex','FontWeight','normal')
ylabel('$y$ ($\mu$m)')
PlotAspect
ax1=nexttile;
C2=[0.87 0.49 0];
C1=[0.95 0.9 0.9];
Cmap=C1+(0:100)'/100.*(C2-C1);
RhoKymo=reshape(AllRho(round(Nx/2),:,:),Nx,size(AllRho,3))';
ActKymo=reshape(AllActin(round(Nx/2),:,:),Nx,size(AllActin,3))';
imagesc((0:Nx-1)*pxlSize,(0:size(RhoKymo,1)-1)*FrTime,RhoKymo)
clim([min(AllRho(:)) max(AllRho(:))])
xlabel('$x$ ($\mu$m)')
ylabel('$t$ (s)')
title('Kymograph','Interpreter','tex','FontWeight','normal')
pbaspect([1 1 1])
colormap(ax1,sky)
ax2=nexttile;
imagesc(ResampledX,ResampledT,DataXCor)
%imagesc((0:Nx-1)*pxlSize,(0:size(RhoKymo,1)-1)*FrTime,ActKymo)
%clim([min(AllActin(:)) max(AllActin(:))])
clim([-1 1])
xlabel('$\Delta r$ ($\mu$m)')
ylabel('$\Delta t$ (s)')
title('$R_{\rho f}$')
%ylabel('$t$ (s)')
pbaspect([1 1 1])
%yticklabels('')
%xticklabels('')
colormap(ax2,turbo)
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
[~,indForT]=min(abs(TimeLags-4));
DataRhoACor = ACorsRho(1:nModesForClassif,1:nModesForClassif,indForT);
DataActACor = ACorsAct(1:nModesForClassif,1:nModesForClassif,indForT);