% This file simulates the orientational actin model in Fig. 5(E,F) in the
% paper. The idea is that actin filaments are transported by advection
% rather than diffusion. We keep track of the density of actin
% f(x,t,theta), where theta is the tangent vector angle in 2D
%function [Statistics,st] = RhoAndActinTauPDEs(Params,seed,postproc)
load('Params.mat')
seed=3;
postproc=1;
rng(seed);
Params=pShortActin; % For Fig. 5E
Params=pActin2; % For Fig. 5F
MakeMovie=1;
advorder=1;
RandomOrient=1;
koff0=Params(1);
rf = Params(2);
Nuc0=Params(3);
NucEn=Params(4);
koffAct = Params(5);
Du = 0.1; % The size of the waves depends on Du
vp = 1;%Params(6);
kbasal=Params(7);
kfb=Params(8);
KFB=Params(9);
ChangeEvery = 1/koffAct;
L=20;
Nx=100; 
Nth = 16;
% Solve for the steady states (no transport or diffusion)
[rts,rtstab,mmg,nmg] = PDERoots(Params,Du,L,Nx);
if (isscalar(rts(:,1)))
    ICRange = [0 2*rts(:,1); 0 2*rts(:,2)];
else 
    ICRange = [min(rts(:,1)) max(rts(:,1)); ...
        min(rts(:,2)) max(rts(:,2))];
end

dt = 0.025; % Stability limit is 1
tf = 240;
tsaves = [40];
saveEvery=floor(1e-6+1/dt);
ChangeEverySteps = ceil(ChangeEvery/dt);

% Initialize grid
dx=L/Nx;
x=(0:Nx-1)*dx;
y=(0:Nx-1)*dx;
[xg,yg]=meshgrid(x,y);
dth = 2*pi/Nth;
th = (1/2:Nth)*dth;
% Set up a wave initial condition (works for worms)
if (Params(6)>Du)
u = ICRange(1,1)+0.5*(1+sin(mmg*2*pi*xg/L).*sin(nmg*2*pi*yg/L))...
   *(ICRange(1,2)-ICRange(1,1));
v0 = ICRange(2,1)+0.5*(1+sin(mmg*2*pi*xg/L).*sin(nmg*2*pi*yg/L))...
   *(ICRange(2,2)-ICRange(2,1)); 
else
load('Initials.mat') %(starfish)
u = uIC;
v0 = sum(vIC,3)*2*pi/8;
end
xvels = zeros(Nx,Nx,Nth);
yvels = zeros(Nx,Nx,Nth);
thForAvg = zeros(Nx,Nx,Nth);
v = zeros(Nx,Nx,Nth);
for j=1:Nth
    v(:,:,j)=v0/(2*pi);
    thForAvg(1:Nx,1:Nx,j)=th(j);
    xvels(1:Nx,1:Nx,j)=vp*cos(th(j));
    yvels(1:Nx,1:Nx,j)=vp*sin(th(j));
end
kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
[kx,ky]=meshgrid(kvals);
ksq=kx.^2+ky.^2;
DivFacFourier_u = (1/dt+Du*ksq);
nSt = floor(tf/dt+1e-10);
nSave = floor(tf/(saveEvery*dt)+1e-10);
AllRho=zeros(Nx,Nx,nSave);
AllActin=zeros(Nx,Nx,nSave);
AllAngles=zeros(Nx,Nx,nSave);

if (MakeMovie)
    close all;
    f=figure('Position',[100 100 700 400]);
end

% Run simulation
PlotUs=zeros(Nx^2,length(tsaves));
PlotVs=zeros(Nx^2,length(tsaves));
PlotTs=0*tsaves;
st=1;
for iT=1:nSt
    vTot = sum(v,3)*dth;
    % Evolve u
    RHS_u = (kbasal+kfb*u.^3./(KFB+u.^3))-(koff0+rf*vTot).*u;
    RHSHat_u = fft2(u/dt+RHS_u);
    uHatNew = RHSHat_u./DivFacFourier_u;
    uNew = ifft2(uHatNew);
    % Evolve v (each theta section separately)
    vNew = zeros(Nx,Nx,Nth);
    % Partition up nucleation rates 
    if (RandomOrient && mod(iT-1,ChangeEverySteps)==0)
        % Each spatial region gets a random theta for nucleation
        Nuc0s = zeros(Nx,Nx,Nth);
        NucEns = zeros(Nx,Nx,Nth);
        % Make rates uniform in theta but not x
        SquareRegionSize = 16; % in um^2
        Nreg = ceil(L^2/SquareRegionSize);
        Regions = MatrixPartition(Nreg,L,Nx,dx);
        Angles0ByRegion = ceil(rand(Nreg,1)*Nth);
        Angles1ByRegion = ceil(rand(Nreg,1)*Nth);
        for iX=1:Nx
            for iY=1:Nx
                Nuc0s(iY,iX,Angles0ByRegion(Regions(iY,iX)))...
                    =Nuc0/dth;
                NucEns(iY,iX,Angles1ByRegion(Regions(iY,iX)))...
                    =NucEn/dth;
            end
        end
    elseif (~RandomOrient)
        Nuc0s = Nuc0*ones(Nx,Nx,Nth)/(2*pi);
        NucEns = NucEn*ones(Nx,Nx,Nth)/(2*pi);
    end
    for j=1:Nth
        RHS_v = (Nuc0s(:,:,j)+NucEns(:,:,j).*u.^2) ...
           - koffAct*v(:,:,j);
        advTerm = advectiveTerm2D(xvels(:,:,j),yvels(:,:,j),...
            v(:,:,j),dx,advorder);
        vNew(:,:,j) = v(:,:,j) + dt*(-advTerm+RHS_v);
    end
    u = uNew;
    v = vNew;
    if (max(abs(u(:)) > 1e5))
        st=0;
        warning('Rejecting because of unstable simulation')
        Statistics.XCor=0;
        Statistics.MeanActin=0;
        return;
    end
    if (mod(iT,saveEvery)==0 && MakeMovie)
        tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
        ax1=nexttile;
        imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
        set(gca,'YDir','Normal')
        %clim([0 StSt(end)])
        %clim([0 0.9])
        %colormap("turbo")
        title(sprintf('Rho; $t= %1.1f$',iT*dt-40))
        %clim(ICRange(1,:))
        colormap(ax1,"sky")
        %colorbar
        pbaspect([1 1 1])
        xlabel("$x$ ($\mu$m)")
        ylabel("$y$ ($\mu$m)")
        ax2=nexttile;
        vBarX = sum(v,3)*dth;
        imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,vBarX);
        set(gca,'YDir','Normal')
        pbaspect([1 1 1])
        C2=[0.87 0.49 0];
        C1=[0.95 0.9 0.9];
        Cmap=C1+(0:100)'/100.*(C2-C1);
        colormap(ax2,Cmap)
        vBarThet=sum(v.*thForAvg*dth,3)./vBarX;
        [~,ThetInd] = max(v,[],3);
        vBarThet=th(ThetInd);
        hold on
        quiver(xg(1:5:end,1:5:end),yg(1:5:end,1:5:end),...
            cos(vBarThet(1:5:end,1:5:end)),...
            sin(vBarThet(1:5:end,1:5:end)),0.5,'k')
        %colorbar
        title(sprintf('Actin; $t= %1.1f$',iT*dt-40))
        xlabel("$x$ ($\mu$m)")
        %clim([6 20])
        %clim(ICRange(2,:))
        movieframes(iT)=getframe(f);
    end
    if (mod(iT-1,saveEvery)==0)
        AllRho(:,:,(iT-1)/saveEvery+1)=u;
        vBarX = sum(v,3)*dth;
        AllActin(:,:,(iT-1)/saveEvery+1)=vBarX;
        vBarThet=sum(v.*thForAvg*dth,3)./vBarX;
        [~,ThetInd] = max(v,[],3);
        vBarThet=th(ThetInd);
        AllAngles(:,:,(iT-1)/saveEvery+1)=vBarThet;
    end
end
%return

    if (postproc)
        AllActin=AllActin(:,:,41:end);
        AllRho=AllRho(:,:,41:end);
        AllAngles=AllAngles(:,:,41:end);
        if (size(rts(:,1))>1)
            Thres=rts(1,2);
        else
            Thres=mean(AllRho(:));
        end
        RhoThres=AllRho>Thres;
        [~,~,nFr]=size(RhoThres);
        NumExcitations=zeros(nFr,1);
        ExSizes=[];
        for iT=1:nFr
            CC = bwconncomp(RhoThres(:,:,iT));
            L2=CC2periodic(CC,[1 1],'L');
            NumExcitations(iT)=max(L2(:));
            for iJ=1:max(L2(:))
                ExSizes=[ExSizes;sum(L2(:)==iJ)/(Nx^2)*L^2];
            end
        end
        [rSim,tSim,XCorsSim] = CrossCorrelations(dx,dx,dt*saveEvery,...
            AllRho,AllActin,0);
        Statistics.XCor=XCorsSim/max(abs(XCorsSim(:)));
        Statistics.rSim=rSim;
        Statistics.tSim=tSim;
        Statistics.ExSizes=ExSizes;
        Statistics.NumExcitations=NumExcitations;
        Statistics.MeanActin=mean(AllActin(:));
    else
        Statistics=[];
    end
%end

% Make plots
tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact')
nexttile
imagesc(x,y,AllRho(:,:,1))
clim([min(AllRho(:)) max(AllRho(:))])
pbaspect([1 1 1])
ylabel('$y$ ($\mu$m)')
set(gca,'YDir','Normal')
nexttile
imagesc(x,y,AllActin(:,:,1))
clim([min(AllActin(:)) max(AllActin(:))])
hold on
quiver(xg(1:5:end,1:5:end),yg(1:5:end,1:5:end),...
cos(AllAngles(1:5:end,1:5:end,1)),...
sin(AllAngles(1:5:end,1:5:end,1)),'w','LineWidth',1)
pbaspect([1 1 1])
set(gca,'YDir','Normal')
yticklabels('')
colormap sky
RhoKymo=reshape(AllRho(75,:,1:end),Nx,size(AllRho,3))'; 
nexttile
imagesc(x,0:size(AllActin,3)-1,RhoKymo)
clim([min(AllRho(:)) max(AllRho(:))])
pbaspect([1 1 1])
ylabel('$t$ (s)')
xlabel('$x$ ($\mu$m)')
ActKymo=reshape(AllActin(75,:,1:end),Nx,size(AllRho,3))';
% Quiver on top of that
u = reshape(cos(AllAngles(75,:,1:end)),Nx,size(AllRho,3))';
v = reshape(sin(AllAngles(75,:,1:end)),Nx,size(AllRho,3))';
[xpl,tpl]=meshgrid(x,(0:size(AllActin,3)-1)/9);
nexttile
imagesc(x,(0:size(AllActin,3)-1)/9,ActKymo)
clim([min(AllActin(:)) max(AllActin(:))])
hold on
quiver(xpl(1:5:end,1:5:end),tpl(1:5:end,1:5:end),...
    u(1:5:end,1:5:end),0*v(1:5:end,1:5:end),'w','LineWidth',1)
pbaspect([1 1 1])
yticklabels('')
xlabel('$x$ ($\mu$m)')
% nexttile
% imagesc(rSim,tSim,Statistics.XCor)
% xlim([0 5])
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