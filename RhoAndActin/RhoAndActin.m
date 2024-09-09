% % In 2D, we are stacking columns
% rng(1);
% kbasal=0.05;
% kfb=1;
% KFB=0.1;
% koff0=0.45;
% rf = 1.25;
% dt=0.1; % Stability limit is 0.4
% tf = 200;
% Du=1e-1; % The size of the waves depends on Du
% tsaves = [15 30 50 100 200];
% % Parameters for the actin
% PoreSize=[];
% ds=0.1;
% FullLifetime = 8; % the average lifetime of a filament (s)
% GrowAmt = floor(1.5/ds*dt);% (#mon per time step - first number is um/s)
% ShrinkAmt = floor(3/ds*dt); % (#mon per time step - first number is um/s)
% MaxLength = 6;
% ICScale = 2.5; 
% gw=0.1;
% ForceFactor=0.4/gw;
% nFilSt = 100;
% ActinUpdate=1;
% saveEvery=10;
% 
% L=20;
% Nx=L/gw; % The grid spacing (should be about gw, always check)
% dx=L/Nx;
% x=(0:Nx-1)*dx;
% y=(0:Nx-1)*dx;
% [xg,yg]=meshgrid(x,y);
% u = rand(Nx,Nx)*ICScale;
% kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
% [kx,ky]=meshgrid(kvals);
% ksq=kx.^2+ky.^2;
% DivFacFourier = (1/dt+Du*ksq);
% nSt = floor(tf/dt+1e-10);
% 
% % If we are applying a stimulus
% nModesStim = Nx;
% Stim = rand(Nx,Nx);
% % Filter out high wave #s
% kmax=(nModesStim*2*pi/L)^2;
% NoiseHat=fft2(Stim);
% NoiseHat(ksq>kmax)=0;
% Stim=ifft2(NoiseHat);
% StimAmp = 0.7;
% Stim=StimAmp*(Stim-min(Stim(:)))./(max(Stim(:))-min(Stim(:)));
% f=figure;
% % Set up an initial actin grid
% [Xf,nPerFil,TimeAtMax] = SetUpActin(L,ds,[],PoreSize,nFilSt,MaxLength);
% if (~isempty(Xf))
%     fg = SpreadToGrid(x,y,Xf,ForceFactor*ds*ones(sum(nPerFil),1),gw); 
% else
%     fg=zeros(Nx);
% end
% 
% % Run simulation
% uv1=zeros(nSt,1);
% PlotUs=zeros(Nx^2,length(tsaves));
% PlotXfs=cell(length(tsaves));
% PlotTs=0*tsaves;
% for iT=1:nSt
%     t=(iT-1)*dt;
%     if (ActinUpdate)
%         % Update old filaments
%         [Xf,nPerFil,TimeAtMax,AddedPts,RemovedPts] = UpdateActin(Xf,nPerFil,...
%             t,TimeAtMax,FullLifetime,GrowAmt,ShrinkAmt,MaxLength,ds);
%         % Nucleate new filaments (proportional to Rho)
%         NucRate=0.002*u.^2;
%         AddedNucs=NucleateNewActin(NucRate,x,y);
%         [nNucs,~]=size(AddedNucs);
%         AddedPts=[AddedPts;AddedNucs];
%         Xf=[Xf;AddedNucs];
%         TimeAtMax=[TimeAtMax;inf*ones(nNucs,1)];
%         nPerFil=[nPerFil;ones(nNucs,1)];
%         % Spread to grid
%         if (~isempty(RemovedPts) || ~isempty(AddedPts))
%             [nRem,~]=size(RemovedPts);
%             [nAdd,~]=size(AddedPts);
%             fgDiff = SpreadToGrid(x,y,[AddedPts;RemovedPts],...
%                 [ForceFactor*ds*ones(nAdd,1); -ForceFactor*ds*ones(nRem,1)],gw);
%             fg=fg+fgDiff;
%         end
%     end
%    % uc = 1-sum(u)*dx^2;
%     RHS = (kbasal+kfb*u.^3./(KFB+u.^3))  - (koff0+rf*fg).*u;
%     RHSHat = fft2(u/dt+RHS);
%     uHatNew = RHSHat./(DivFacFourier);
%     uNew = ifft2(uHatNew);
%     u = uNew;
%     if (sum(abs(iT*dt-tsaves)<1e-10)>0)
%         index = find(abs(iT*dt-tsaves)<1e-10);
%         PlotUs(:,index)=reshape(u,[],1);
%         PlotXfs{index}=Xf-floor(Xf/L)*L;
%         PlotTs(index)=iT*dt;
%     end
%     if (mod(iT-1,20)==0)
%         imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
%         set(gca,'YDir','Normal')
%         hold on
%         if (~isempty(Xf))
%             Xfpl=Xf-floor(Xf/L)*L;
%             plot(Xfpl(:,1),Xfpl(:,2),'ko','MarkerSize',0.25)
%         end
%         colorbar
%         clim([0 2.5])
%         colormap("turbo")
%         title(sprintf('$t= %1.1f$',t))
%         drawnow
%         hold off
%         movieframes(iT)=getframe(f);
%     end
%     uv1(iT)=u(Nx/2,Nx/2);
%     TotActin(iT)=sum(nPerFil);
%     if (mod(iT-1,saveEvery)==0)
%         AllRho(:,:,(iT-1)/saveEvery+1)=u;
%         AllActin(:,:,(iT-1)/saveEvery+1)=fg;
%     end
% end
% 
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
%         plot(PlotXfs{iT}(:,1),PlotXfs{iT}(:,2),'ko','MarkerSize',0.2)
%     end
%     if (iT==nPlot)
%         colorbar
%     end
%     if (iT==1)
%         ylabel('$y$')
%     end
% end

figure;
padxy=0;
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
% nexttile
% plot((dt:dt:tf),TotActin)
% xlabel('$t$ (s)')
% ylabel('Polymerized actin ($\mu$m)')
% title('Total actin')
nexttile
[Uvals,dtvals,DistsByR] = CrossCorrelations(dx,dx,dt*saveEvery,...
    AllRho,AllRho,padxy);
imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
title('Rho')
ylim([-60 60])
xlim([0 10])
a=clim;
clim([-max(abs(a)) max(abs(a))]);
xlabel('$\Delta r$')
ylabel('$\Delta t$')
nexttile
[Uvals,dtvals,DistsByR] = CrossCorrelations(dx,dx,dt*saveEvery,...
    AllRho,AllActin,padxy);
imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
xlabel('$\Delta r$')
ylim([-60 60])
xlim([0 10])
title('Rho-Actin')
a=clim;
clim([-max(abs(a)) max(abs(a))]);
nexttile
[Uvals,dtvals,DistsByR] = CrossCorrelations(dx,dx,dt*saveEvery,...
    AllActin,AllActin,padxy);
imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
xlabel('$\Delta r$')
title('Actin')
ylim([-60 60])
xlim([0 10])
a=clim;
clim([-max(abs(a)) max(abs(a))]);
colormap turbo
