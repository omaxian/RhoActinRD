nSamp = 1000;
k0s=0.9;
rfs=0.4;
FullLifetime=15;
GrShrnk = 0.1; % rate in um/s
MaxLength = 1; % in um
N0PerMax = 0.2;
Du = 1e-1;
Params = [k0s;rfs;FullLifetime;GrShrnk;MaxLength;N0PerMax;Du];
DifferencesModelExp=zeros(nSamp,1);
AllParameters=zeros(7,nSamp);
StableSim = ones(nSamp,1);
LastAccept=1;
Accepted=zeros(nSamp,1);
%inds=23;
for iSamp=1:nSamp
% Params=AllParameters(:,inds(iSamp));
% DiffChk=DifferencesModelExp(inds(iSamp))
OldParams=Params;
if (iSamp>1)
    % Proposal
    k0Step = min(1.25-OldParams(1),0.05*randn);
    rfStep = 0.05*randn;
    LifetimeStep =2*randn;
    GrShrnkStep = 0.05*randn;
    MaxLenStep = 0.2*randn;
    N0PerMaxStep=0.04*randn; % symmetric proposal
    DuStep = 0.01*randn;
    DeltaParams=max([k0Step;rfStep;LifetimeStep;...
        GrShrnkStep;MaxLenStep;N0PerMaxStep;DuStep],-Params);
    Params = OldParams+DeltaParams;
end
close all;
rng(2);
kbasal=0.05;
kfb=1;
KFB=0.1;
kNoise=0;
koff0=Params(1);
rf = Params(2);
% Solve for the steady states in absence of actin
try
p=0:0.001:100;
OnRate = (kbasal+kfb*p.^3./(KFB+p.^3));
OffRate=(koff0*p);
Net=OnRate-OffRate;
SgnChg=find((Net(1:end-1).*Net(2:end))<0);
StSt=p(SgnChg);
Nuc0=Params(6)/max(StSt);
catch
StSt=100;
Nuc0=Params(6)/max(StSt);
end
dt=0.1; % Stability limit is 0.4
tf = 240;
Du=Params(7); % The size of the waves depends on Du
tsaves = [100 200];
% Parameters for the actin
PoreSize=[];
ds=0.1;
FullLifetime = Params(3); % the average lifetime of a filament (s)
GrShrnk=Params(4);
GrowAmt = (GrShrnk/ds*dt);% (#mon per time step - first number is um/s)
ShrinkAmt = (GrShrnk/ds*dt); % (#mon per time step - first number is um/s)
MaxLength = max(ds,floor(Params(5)/ds)*ds);
ICScale = StSt(end); 
gw=0.1;
ForceFactor=0.4/gw;
nFilSt = ceil(500/MaxLength);
ActinUpdate=1;
saveEvery=10;

L=20;
Nx=L/gw; % The grid spacing (should be about gw, always check)
dx=L/Nx;
x=(0:Nx-1)*dx;
y=(0:Nx-1)*dx;
[xg,yg]=meshgrid(x,y);
u = rand(Nx,Nx)*ICScale;
kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
[kx,ky]=meshgrid(kvals);
ksq=kx.^2+ky.^2;
DivFacFourier = (1/dt+Du*ksq);
nSt = floor(tf/dt+1e-10);

% If we are applying a stimulus
nModesStim = Nx;
Stim = rand(Nx,Nx);
% Filter out high wave #s
kmax=(nModesStim*2*pi/L)^2;
NoiseHat=fft2(Stim);
NoiseHat(ksq>kmax)=0;
Stim=ifft2(NoiseHat);
StimAmp = 0.7;
Stim=StimAmp*(Stim-min(Stim(:)))./(max(Stim(:))-min(Stim(:)));
f=figure;
% Set up an initial actin grid
[Xf,nPerFil,~] = SetUpActin(L,ds,[],PoreSize,nFilSt,MaxLength);
TimeAtMax = -1.0*rand(length(nPerFil),1)*FullLifetime;
if (~isempty(Xf))
    fg = SpreadToGrid(x,y,Xf,ForceFactor*ds*ones(sum(nPerFil),1),gw); 
else
    fg=zeros(Nx);
end

% Run simulation
uv1=zeros(nSt,1);
PlotUs=zeros(Nx^2,length(tsaves));
PlotXfs=cell(length(tsaves));
PlotTs=0*tsaves;
for iT=1:nSt
    t=(iT-1)*dt;
    if (ActinUpdate)
        % Update old filaments
        [Xf,nPerFil,TimeAtMax,AddedPts,RemovedPts] = UpdateActin(Xf,nPerFil,...
            t,TimeAtMax,FullLifetime,GrowAmt,ShrinkAmt,MaxLength,ds);
        % Nucleate new filaments (proportional to Rho)
        NucRate=Nuc0*u.^2*dt;
        AddedNucs=NucleateNewActin(NucRate,x,y);
        [nNucs,~]=size(AddedNucs);
        AddedPts=[AddedPts;AddedNucs];
        Xf=[Xf;AddedNucs];
        TimeAtMax=[TimeAtMax;inf*ones(nNucs,1)];
        nPerFil=[nPerFil;ones(nNucs,1)];
        % Spread to grid
        if (~isempty(RemovedPts) || ~isempty(AddedPts))
            [nRem,~]=size(RemovedPts);
            [nAdd,~]=size(AddedPts);
            fgDiff = SpreadToGrid(x,y,[AddedPts;RemovedPts],...
                [ForceFactor*ds*ones(nAdd,1); -ForceFactor*ds*ones(nRem,1)],gw);
            fg=fg+fgDiff;
        end
    end
    RHS = (kbasal+kfb*u.^3./(KFB+u.^3))-(koff0+rf*fg).*u;
    RHSHat = fft2(u/dt+RHS);
    uHatNew = RHSHat./(DivFacFourier);
    uNew = ifft2(uHatNew);
    u = uNew;
    if (max(abs(u(:)) > 1e5))
        StableSim(iSamp)=0;
        break
    end
    if (sum(abs(iT*dt-tsaves)<1e-10)>0)
        index = find(abs(iT*dt-tsaves)<1e-10);
        PlotUs(:,index)=reshape(u,[],1);
        PlotXfs{index}=Xf-floor(Xf/L)*L;
        PlotTs(index)=iT*dt;
    end
    if (mod(iT-1,5)==0)
        imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
        set(gca,'YDir','Normal')
        hold on
        if (~isempty(Xf))
            Xfpl=Xf-floor(Xf/L)*L;
            plot(Xfpl(:,1),Xfpl(:,2),'ko','MarkerSize',0.25)
        end
        colorbar
        clim([0 StSt(end)])
        colormap("turbo")
        title(sprintf('$t= %1.1f$',t))
        drawnow
        hold off
        movieframes(iT)=getframe(f);
    end
    if (mod(iT-1,saveEvery)==0)
        AllRho(:,:,(iT-1)/saveEvery+1)=u;
        AllActin(:,:,(iT-1)/saveEvery+1)=fg;
    end
end
try
CompareXCors;
catch
DiffNorm=inf;
InterpolatedSim=[];
end
DifferencesModelExp(iSamp)=DiffNorm;
AllParameters(:,iSamp)=Params;
AllXCors{iSamp}=InterpolatedSim;
% Accept or reject
if (iSamp > 1)
    kT=3.7;
    r1=exp(-DifferencesModelExp(iSamp)/kT);
    rLast=exp(-DifferencesModelExp(LastAccept)/kT);
    pAcc=r1/rLast;
    if (rand > pAcc) % reject
        Params = OldParams;
    else
        LastAccept=iSamp;
        Accepted(iSamp)=1;
    end
end
end
% figure;
% [~,nPlot]=size(PlotUs);
% tiledlayout(1,nPlot+1,'Padding', 'none', 'TileSpacing', 'compact');
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
% nexttile
% imagesc(rSim,tSim,XCorsSim/max(abs(XCorsSim(:))))
% clim([-1 1])
% colorbar
% % 
% figure;
% padxy=0;
% tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
% % nexttile
% % plot((dt:dt:tf),TotActin)
% % xlabel('$t$ (s)')
% % ylabel('Polymerized actin ($\mu$m)')
% % title('Total actin')
% nexttile
% [Uvals,dtvals,DistsByR] = CrossCorrelations(dx,dx,dt*saveEvery,...
%     AllRho,AllRho,padxy);
% imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
% title('Rho')
% ylim([-60 60])
% xlim([0 10])
% a=clim;
% clim([-max(abs(a)) max(abs(a))]);
% xlabel('$\Delta r$')
% ylabel('$\Delta t$')
% nexttile
% [Uvals,dtvals,DistsByR] = CrossCorrelations(dx,dx,dt*saveEvery,...
%     AllRho,AllActin,padxy);
% imagesc(Uvals,dtvals,DistsByR/max(abs(DistsByR(:))))
% xlabel('$\Delta r$')
% ylim([-60 60])
% xlim([0 10])
% title('Rho-Actin')
% a=clim;
% clim([-max(abs(a)) max(abs(a))]);
% nexttile
% [Uvals,dtvals,DistsByR] = CrossCorrelations(dx,dx,dt*saveEvery,...
%     AllActin,AllActin,padxy);
% imagesc(Uvals,dtvals,DistsByR/max(DistsByR(:)))
% xlabel('$\Delta r$')
% title('Actin')
% ylim([-60 60])
% xlim([0 10])
% a=clim;
% clim([-max(abs(a)) max(abs(a))]);
% colormap turbo
