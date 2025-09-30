%figure;
% %tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact')
% ax1=nexttile
% imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,AllRho(:,:,1))
% set(gca,'YDir','Normal')
% pbaspect([1 1 1])
% ylabel('$y$ ($\mu$m)')
% xticklabels('')
% clim([min(AllRho(:)) max(AllRho(:))])
% colormap(ax1,sky)
% ax3=nexttile;
% imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,AllActin(:,:,1))
% set(gca,'YDir','Normal')
% pbaspect([1 1 1])
% C2=[0.87 0.49 0];
% C1=[0.95 0.9 0.9];
% Cmap=C1+(0:100)'/100.*(C2-C1);
% colormap(ax3,Cmap)
% %ylabel('$y$ ($\mu$m)')
% xticklabels('')
% clim([min(AllActin(:)) max(AllActin(:))])
% yticklabels('')
ax4=nexttile;
RhoKymo=reshape(AllRho(3*Nx/4,:,1:end),Nx,size(AllRho,3))';
ActKymo=reshape(AllActin(3*Nx/4,:,1:end),Nx,size(AllActin,3))';
imagesc((0:Nx-1)*dx,0:size(RhoKymo,3)-1,RhoKymo)
clim([min(min(AllRho(:)),rts(1,1)) max(max(AllRho(:)),rts(3,1))])
xlabel('$x$ ($\mu$m)')
ylabel('$t$ (s)')
pbaspect([1 1 1])
colormap(ax4,sky)
title(strcat('$D_f =\,$',num2str(Params(6)),...
    '; $k_\textrm{diss} = \,$',num2str(Params(5))))
% ax2=nexttile;
% imagesc((0:Nx-1)*dx,0:size(RhoKymo,3)-1,ActKymo)
% clim([min(AllActin(:)) max(AllActin(:))])
% xlabel('$x$ ($\mu$m)')
% %ylabel('$t$ (s)')
% pbaspect([1 1 1])
% yticklabels('')
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
%colormap(ax2,Cmap)
