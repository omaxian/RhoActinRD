% Cross correlation profiles in experimental data
MovieNum=3;
pxlSize = 0.1;
FrTime = 0.3;
Rho=double(load(strcat('NMYRho_',num2str(MovieNum),'.mat')).RhoData);
Actin=double(load(strcat('NMYActin_',num2str(MovieNum),'.mat')).ActinData);
% Bement data
% pxlSize = 120/212;
% FrTime = 330/72;
% Rho=tiffreadVolume("BementRho.tif");
% Actin=tiffreadVolume("BementRGA.tif");
% %Difference subtraction (removes static signal)
% Rho=Rho-Rho(:,:,end);
% Actin=Actin-Actin(:,:,end);
% Filter data
nModes=inf;
Rho = FilterData(Rho,nModes);
Actin = FilterData(Actin,nModes);
% Make a plot of the data
%tsPl=[0:4:20];
tsPl=[10 20 30];
FrPl=ceil(tsPl/FrTime)+1;
%Nx=length(xc);
Nx=200;
x=(0:Nx-1)*pxlSize;
y=(0:Nx-1)*pxlSize;
tiledlayout(2,length(tsPl),'Padding', 'none', 'TileSpacing', 'compact');
for iP=1:length(tsPl)
    nexttile
    imagesc(x,y,Rho(:,:,FrPl(iP)));
    title(sprintf('$t= %1.1f$',FrTime*(FrPl(iP)-1)))
    if (iP==1)
        ylabel('Rho concentration')
    end
end
for iP=1:length(tsPl)
    nexttile
    imagesc(x,y,Actin(:,:,FrPl(iP)));
    if (iP==1)
        ylabel('Actin concentration')
    end
end
colormap turbo


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