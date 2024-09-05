%close all;
% In 2D, we are stacking columns
rng(1);
L=1;
kbasal=0.05;
kfb=1;
KFB=0.1;
koff0=1;
rf = 0.5;
dt=0.2; % Stability limit is 0.4
tf = 25;
Nx=200; % The grid spacing (always check)
Du=2.5e-4; % The size of the waves depends on Du
tsaves = [50 100 150 225 300];
ICScale = 0.1;

dx=L/Nx;
x=(0:Nx-1)*dx;
y=(0:Nx-1)*dx;
[xg,yg]=meshgrid(x,y);
kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
[kx,ky]=meshgrid(kvals);
ksq=kx.^2+ky.^2;
DivFacFourier = (1/dt+Du*ksq);

% Filter out high wave #s
nModesFilt = 10;
kFilt=(nModesFilt*2*pi/L)^2;

% Actin set up (probably should be periodically extending this data)
MovieNum=3;
StFrame=230;
FrTime = 0.3;
Rho=double(load(strcat('NMYRho_',num2str(MovieNum),'.mat')).RhoData);
Actin=double(load(strcat('NMYActin_',num2str(MovieNum),'.mat')).ActinData);
[~,~,nFr]=size(Actin);
%tf=(nFr-1)*dt;
Actin(:,:,1:StFrame-1)=[];
Rho(:,:,1:StFrame-1)=[];

% Filter the data
FilteredActin = FilterData(Actin,ksq,kFilt,Nx);
FilteredRho = FilterData(Rho,ksq,kFilt,Nx);
% Get everything on [0,1]
FilteredActin = FilteredActin/max(FilteredActin(:));
FilteredRho = FilteredRho/max(FilteredRho(:));

u = 0*FilteredRho(:,:,1);
nSt = floor(tf/dt+1e-10);
uv1=zeros(nSt,1);
f=figure('Position',[100 200 1300 600]);
for iT=1:nSt
    xc=45:65;
    yc=55:75;
    t=(iT-1)*dt;
    ExpIndex = 1+floor((iT-1)*dt/FrTime);
    uv1(iT,1)=mean(mean(FilteredRho(yc,xc,ExpIndex)));
    uv1(iT,2)=mean(mean(FilteredActin(yc,xc,ExpIndex)));
    uv1(iT,3)=mean(mean(u(yc,xc)));

    fg=FilteredActin(:,:,ExpIndex); % Filter

    RHS = (kbasal+kfb*u.^3./(KFB+u.^3))  - (koff0+rf.*(fg>0.85)).*u;
    % Stimulate the box
    if (t > 7 && t < 8.5)
        RHS(yc,xc)=RHS(yc,xc)+1;
    end
    RHSHat = fft2(u/dt+RHS);
    uHatNew = RHSHat./(DivFacFourier);
    uNew = ifft2(uHatNew);
    u = uNew;
    if (mod(iT-1,4)==0)
        iSave = (iT-1)/4+1;
        tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
        nexttile
        imagesc(FilteredRho(1:Nx/2,1:Nx/2,ExpIndex))
        clim([min(FilteredRho(:)) max(FilteredRho(:))])
        PlotAspect
        colorbar
        hold on
        rectangle('Position',[min(xc) min(yc) range(xc) range(yc)])
        title({sprintf('$t=%1.1f$, Fr=%d',t,ExpIndex)},{'Experimental Rho'})
        nexttile
        imagesc(u(1:Nx/2,1:Nx/2))
        clim([0 1])
        colorbar
        PlotAspect
        hold on
        rectangle('Position',[min(xc) min(yc) range(xc) range(yc)])
        title({sprintf('$t=%1.1f$, Fr=%d',t,ExpIndex)},{'Simulated Rho'})
        nexttile
        imagesc(FilteredActin(1:Nx/2,1:Nx/2,ExpIndex))
        clim([min(FilteredActin(:)) max(FilteredActin(:))])
        title({sprintf('$t=%1.1f$, Fr=%d',t,ExpIndex)},{'Experimental actin'})
        hold on
        rectangle('Position',[min(xc) min(yc) range(xc) range(yc)])
        PlotAspect
        colorbar
        colormap turbo
        movieframes(iSave)=getframe(f);
    end
end

% figure;
% [~,nPlot]=size(PlotUs);
% tiledlayout(1,nPlot,'Padding', 'none', 'TileSpacing', 'compact');
% for iT=1:nPlot
%     nexttile
%     imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,reshape(PlotUs(:,iT),Nx,Nx));
%     title(strcat('$t=$',num2str(PlotTs(iT))))
%     clim([min(PlotUs(:)) max(PlotUs(:))])
%     %clim([0 3])
%     set(gca,'YDir','Normal')
%     colormap(turbo)
%     hold on
%     if (~isempty(Xf))
%         plot(PlotXfs(1:length(Xf),iT),PlotXfs(length(Xf)+1:end,iT),'ko','MarkerSize',0.2)
%     end
%     if (iT==nPlot)
%         colorbar
%     end
%     if (iT==1)
%         ylabel('$y$')
%     end
% end
figure
plot(0:dt:tf-dt,uv1)
title('Dynamics in the box')