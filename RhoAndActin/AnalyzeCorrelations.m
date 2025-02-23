% Cross correlation profiles in experimental data
Name='nmy-cyk';
Widths = [50];
Threses = [0.1];
TimePlot=400;
for iFw=1:length(Widths)
%nexttile
for iThr=1:length(Threses)
MovieNum=3;
pxlSize = 0.1;
Rho=double(load(strcat(Name,'Rho_',num2str(MovieNum),'.mat')).RhoData);
Actin=double(load(strcat(Name,'Actin_',num2str(MovieNum),'.mat')).ActinData);
Info=load(strcat(Name,'Info_',num2str(MovieNum),'.mat'));
FName=Info.Info.Name;
FrTime = Info.Info.TimeInt;
% Bement data
% pxlSize = 120/212;
% FrTime = 330/72;
% Rho=tiffreadVolume("BementRho.tif");
% Actin=tiffreadVolume("BementRGA.tif");
% %Difference subtraction (removes static signal)
% Rho=Rho-Rho(:,:,end);
% Actin=Actin-Actin(:,:,end);
% Adjust so that each frame has same mean
GlobalMeanRho=mean(Rho(:));
GlobalMeanActin=mean(Actin(:));
[ny,nx,nFr]=size(Rho);
for iT=1:nFr
    Rho(:,:,iT)=Rho(:,:,iT)-mean(mean(Rho(:,:,iT)))+GlobalMeanRho;
    Actin(:,:,iT)=Actin(:,:,iT)-mean(mean(Actin(:,:,iT)))+GlobalMeanActin;
end
Nx=200;
x=(0:Nx-1)*pxlSize;
y=(0:Nx-1)*pxlSize;
nexttile
imagesc(x,y,Rho(:,:,TimePlot))
xlabel('$x$ ($\mu$m)')
ylabel('$y$ ($\mu$m)')
colormap turbo
% Filter data
%nModes=10;
%nModesTime=10;
%FiltRho = FilterData(Rho,Widths(iFw),Widths(iFw));
Filtx=smoothdata(Rho,1,'sgolay',Widths(iFw));
Filtxy=smoothdata(Filtx,2,'sgolay',Widths(iFw));
FiltRho=smoothdata(Filtxy,3,'sgolay',Widths(iFw));
nexttile
imagesc(x,y,FiltRho(:,:,TimePlot))
xlabel('$x$ ($\mu$m)')
% Identify pulsing regions
Thres=mean(FiltRho(:))+Threses(iThr)*(max(FiltRho(:))-mean(FiltRho(:)));
Excited=FiltRho>Thres;
% tsPl=[20 40 60];
% FrPl=ceil(tsPl/FrTime)+1;

% figure;
% tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
% for iP=1:length(tsPl)
% nexttile
% L2 = bwlabel(Excited(:,:,FrPl(iP)));
% imagesc(x,y,L2);
% title(sprintf('$t= %1.1f$',FrTime*(FrPl(iP)-1)))
% %clim([min(FiltRho(:)) max(FiltRho(:))])
% pbaspect([1 1 1])
% colormap turbo
% end
% end
% end
%FiltActin = FilterData(Actin,nModes,nModesTime);
Filtx=smoothdata(Actin,1,'sgolay',Widths(iFw));
Filtxy=smoothdata(Filtx,2,'sgolay',Widths(iFw));
FiltActin=smoothdata(Filtxy,3,'sgolay',Widths(iFw));
NumExes=zeros(nFr,1);
AvgExSize=zeros(nFr,1);
AllExes=[];
% Compute their number and average area
for iT=1:nFr
    L2 = bwlabel(Excited(:,:,iT));
    nEx = max(L2(:));
    ExSize=zeros(nEx,1);
    for j=1:nEx
        ExSize(j)=sum(sum(L2==j))*pxlSize^2;
    end
    AllExes=[AllExes;ExSize];
    NumExes(iT)=nEx;
    AvgExSize(iT)=mean(ExSize);
end
nexttile
imagesc(x,y,bwlabel(Excited(:,:,TimePlot)))
xlabel('$x$ ($\mu$m)')% Cross correlation
[UvalsF,dtvalsF,DistsByRF] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    FiltRho,FiltActin,1);
DistsByRF=DistsByRF/max(DistsByRF(:));
DistsByRF=DistsByRF(abs(dtvalsF)<tmax,abs(UvalsF)<rmax)/max(DistsByRF(:));
dtvalsF=dtvalsF(abs(dtvalsF)<tmax);
UvalsF=UvalsF(abs(UvalsF)<rmax);
WtsByR = exp(-UvalsF'/2);
WtsByT = exp(-abs(dtvalsF)'/60);
TotWts=WtsByR.*WtsByT;
XCorNorm=TotWts.*DistsByRF.^2;
XCorNorm=sum(XCorNorm(:))
end
end
% figure;
% imagesc(UvalsF,dtvalsF,DistsByRF/max(DistsByRF(:)))
% xlabel('$\Delta r$')
% ylim([-tmax tmax])
% xlim([0 rmax])
% title('Rho-Actin')% x=(0:nx-1)*pxlSize;
% y=(0:ny-1)*pxlSize;
% for iT=1:10:nFr
% imagesc(x,y,FiltRho(:,:,iT))
% hold on
% [xg,yg]=meshgrid(x,y);
% ex=1.0*Excited(:,:,iT);
% ex(ex==0)=1e-16;
% scatter(xg(:),yg(:),1*ex(:),'filled');
% clim([min(FiltRho(:)) max(FiltRho(:))])
% drawnow;
% hold off
% end
% Make a plot of the data
%tsPl=[0:4:20];
return

figure;
padxy=1;
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
[Uvals,dtvals,DistsByR] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    Rho,Rho,padxy);
imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
title('Rho')
ylim([-60 60])
xlim([0 10])
a=clim;
clim([-max(abs(a)) max(abs(a))]);
xlabel('$\Delta r$')
ylabel('$\Delta t$')
nexttile
[Uvals,dtvals,DistsByR] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    Rho,Actin,padxy);
imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
xlabel('$\Delta r$')
ylim([-60 60])
xlim([0 10])
title('Rho-Actin')
a=clim;
clim([-max(abs(a)) max(abs(a))]);
nexttile
[Uvals,dtvals,DistsByR] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    Actin,Actin,padxy);
imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
xlabel('$\Delta r$')
title('Actin')
ylim([-60 60])
xlim([0 10])
a=clim;
clim([-max(abs(a)) max(abs(a))]);
colormap turbo


figure;
padxy=1;
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
[Uvals,dtvals,DistsByR] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    FiltRho,FiltRho,padxy);
imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
title('FRho')
ylim([-60 60])
xlim([0 10])
a=clim;
clim([-max(abs(a)) max(abs(a))]);
xlabel('$\Delta r$')
ylabel('$\Delta t$')
nexttile
[UvalsF,dtvalsF,DistsByR_Filtered] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    FiltRho,FiltActin,padxy);
imagesc(UvalsF,dtvalsF,DistsByR_Filtered/max(DistsByR_Filtered(:)))
xlabel('$\Delta r$')
ylim([-60 60])
xlim([0 10])
title('Rho-Actin')
a=clim;
clim([-max(abs(a)) max(abs(a))]);
nexttile
[Uvals,dtvals,DistsByR] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
    FiltActin,FiltActin,padxy);
imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
xlabel('$\Delta r$')
title('Actin')
ylim([-60 60])
xlim([0 10])
a=clim;
clim([-max(abs(a)) max(abs(a))]);
colormap turbo


%% This stuff is all just for testing
% 1D cross correlation (periodic)
% x=-1:0.001:0.999;
% N=length(x);
% a=exp(-(x/0.1).^2);
% b=circshift(a,-100);
% c=xcorr(a,b);
% c2=circshift(ifft(fft(a).*conj(fft(b))),N/2);
% plot(-N+1:N-1,c)
% hold on
% plot(-N/2:N/2-1,c2)

% 1D non-periodic
% x=0:0.01:1;
% a=exp(-((x-0.3)/0.1).^2);
% b=exp(-((x-0.5)/0.1).^2);
% c=xcorr(a,b);
% N=length(x);
% az=[a zeros(N-1,1)'];
% bz=[b zeros(N-1,1)'];
% c2=circshift(ifft(fft(az).*conj(fft(bz))),N-1);
% plot(-N+1:N-1,c)
% hold on
% plot(-N+1:N-1,c2)

% 2D cross correlation
% x=-1:0.01:0.99;
% Nx=length(x);
% y=-1:0.01:0.99;
% Ny=length(y);
% [xg,yg]=meshgrid(x,y);
% a=exp(-((xg-0.1)/0.1).^2).*exp(-((yg-0.3)/0.1).^2);
% b=exp(-((xg-0.3)/0.1).^2).*exp(-((yg-0.1)/0.1).^2);
% c3=circshift(ifft2(fft2(b).*conj(fft2(a))),[Ny/2 Nx/2]);
% imagesc(x,y,c3)
% 