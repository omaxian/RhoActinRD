tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact')
nexttile
imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,AllRho(:,:,1))
set(gca,'YDir','Normal')
pbaspect([1 1 1])
ylabel('$y$ ($\mu$m)')
xticklabels('')
clim([min(AllRho(:)) max(AllRho(:))])
colormap sky
nexttile
imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,AllActin(:,:,1))
set(gca,'YDir','Normal')
pbaspect([1 1 1])
%ylabel('$y$ ($\mu$m)')
xticklabels('')
clim([min(AllActin(:)) max(AllActin(:))])
yticklabels('')
nexttile
RhoKymo=reshape(AllRho(3*Nx/4,:,1:end),Nx,tf-140)';
ActKymo=reshape(AllActin(3*Nx/4,:,1:end),Nx,tf-140)';
imagesc((0:Nx-1)*dx,0:size(RhoKymo,3)-1,RhoKymo)
clim([min(AllRho(:)) max(AllRho(:))])
xlabel('$x$ ($\mu$m)')
ylabel('$t$ (s)')
pbaspect([1 1 1])
nexttile
imagesc((0:Nx-1)*dx,0:size(RhoKymo,3)-1,ActKymo)
clim([min(AllActin(:)) max(AllActin(:))])
xlabel('$x$ ($\mu$m)')
%ylabel('$t$ (s)')
pbaspect([1 1 1])
yticklabels('')
% nexttile
% imagesc(Statistics.rSim,Statistics.tSim,Statistics.XCor)
% xlim([0 10])
% ylim([-120 120])
% clim([-1 1])
% xlabel('$\Delta r$ ($\mu$m)')
% ylabel('$\Delta t$ (s)')
% pbaspect([1 1 1])
% nexttile
% dsHist=4;
% xp=histcounts(ExSizes,0:dsHist:400);
% xp=xp/(sum(xp)*dsHist);
% plot(dsHist/2:dsHist:400,xp)
% xlim([0 200])
% xlabel('Excitation size ($\mu$m$^2$)')
% ylabel('pdf')
% pbaspect([1 1 1])
C2=[0.87 0.49 0];
C1=[0.95 0.9 0.9];
Cmap=C1+(0:100)'/100.*(C2-C1);
colormap(Cmap)