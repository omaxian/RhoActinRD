function Statistics = RhoAndActin(Params,seed)
    % Parameters:
    % koff0, rf, FullLifetime, GrShrnk, MaxLength, Nuc0 
    % in that order
    % Output is the difference in the cross correlations compared to
    % experimental data
    rng(seed);
    MakeMovie=0;
    kbasal=0.05;
    kfb=1;
    KFB=0.1;
    koff0=Params(1);
    rf = Params(2);
    % Solve for the steady states in absence of actin
    p=0:0.001:100;
    OnRate = (kbasal+kfb*p.^3./(KFB+p.^3));
    OffRate=(koff0*p);
    Net=OnRate-OffRate;
    SgnChg=find((Net(1:end-1).*Net(2:end))<0);
    StSt=p(SgnChg);
    Nuc0=Params(6)/max(StSt)^2;
    dt = 0.25; % Stability limit is 1
    tf = 200;
    Du=0.1; % The size of the waves depends on Du
    tsaves = [120 140];
    % Parameters for the actin
    PoreSize=[];
    ds=0.1;
    FullLifetime = Params(3); % the average lifetime of a filament (s)
    GrShrnk=Params(4);
    GrowAmt = (GrShrnk/ds*dt);% (#mon per time step - first number is um/s)
    ShrinkAmt = (GrShrnk/ds*dt); % (#mon per time step - first number is um/s)
    MaxLength = max(ds,floor(Params(5)/ds)*ds);
    ICScale = StSt(end); 
    gw=ds;
    ForceFactor=0.4/gw;
    nFilSt = ceil(500/MaxLength);
    ActinUpdate=1;
    saveEvery=floor(1e-6+1/dt);
    
    L=20;
    Nx=L/gw*1/2; % The grid spacing 
    dx=L/Nx;
    x=(0:Nx-1)*dx;
    y=(0:Nx-1)*dx;
    %[xg,yg]=meshgrid(x,y);
    u = ones(Nx,Nx)*ICScale;
    kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
    [kx,ky]=meshgrid(kvals);
    ksq=kx.^2+ky.^2;
    DivFacFourier = (1/dt+Du*ksq);
    nSt = floor(tf/dt+1e-10);
    nSave = floor(tf/(saveEvery*dt)+1e-10);
    AllRho=zeros(Nx,Nx,nSave)+min(StSt);
    AllActin=zeros(Nx,Nx,nSave);

    if (MakeMovie)
        close all;
        f=figure;
    end
    % Set up an initial actin grid
    [Xf,nPerFil,~] = SetUpActin(L,ds,[],PoreSize,nFilSt,MaxLength);
    TimeAtMax = -1.0*rand(length(nPerFil),1)*FullLifetime;
    if (~isempty(Xf))
        fg = SpreadToGrid(x,y,Xf,ForceFactor*ds*ones(sum(nPerFil),1),gw); 
    else
        fg=zeros(Nx);
    end
    
    % Run simulation
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
        % Break out of simulations that aren't doing anything
        if (max(u(:)) < 1.05*min(StSt) && min(StSt) < 0.5)
            Statistics.XCor=0;
            Statistics.MeanActin=0;
            return;
        end
        RHS = (kbasal+kfb*u.^3./(KFB+u.^3))-(koff0+rf*fg).*u;
        RHSHat = fft2(u/dt+RHS);
        uHatNew = RHSHat./(DivFacFourier);
        uNew = ifft2(uHatNew);
        u = uNew;
        if (max(abs(u(:)) > 1e5))
            warning('Rejecting because of unstable simulation')
            Statistics.XCor=0;
	    Statistics.MeanActin=0;
            return;
        end
        if (sum(abs(iT*dt-tsaves)<1e-10)>0)
            index = find(abs(iT*dt-tsaves)<1e-10);
            PlotUs(:,index)=reshape(u,[],1);
            PlotXfs{index}=Xf-floor(Xf/L)*L;
            PlotTs(index)=iT*dt;
        end
        if (mod(iT-1,10)==0 && MakeMovie)
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
            set(gca,'YDir','Normal')
            hold on
            if (~isempty(Xf))
                Xfpl=Xf-floor(Xf/L)*L;
                plot(Xfpl(:,1),Xfpl(:,2),'ko','MarkerSize',0.125)
            end
            colorbar
            %clim([0 StSt(end)])
            %colormap("turbo")
            title(sprintf('$t= %1.1f$',t))
            clim([0 max(StSt)])
            drawnow
            hold off
            movieframes(iT)=getframe(f);
        end
        if (mod(iT-1,saveEvery)==0)
            AllRho(:,:,(iT-1)/saveEvery+1)=u;
            AllActin(:,:,(iT-1)/saveEvery+1)=fg;
        end
    end
    % Post-process to get cross correlations and excitation sizes
    % Compute cross correlation function
    AllActin=AllActin(:,:,41:end);
    AllRho=AllRho(:,:,41:end);
    if (size(StSt)>1)
        Thres=StSt(2);
    else
        Thres=0.5;
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
    Statistics.XCor=XCorsSim/max(XCorsSim(:));
    Statistics.rSim=rSim;
    Statistics.tSim=tSim;
    Statistics.ExSizes=ExSizes;
    Statistics.NumExcitations=NumExcitations;
    Statistics.MeanActin=mean(AllActin(:));
    if (0)
        figure;
        [~,nPlot]=size(PlotUs);
        tiledlayout(1,nPlot+1,'Padding', 'none', 'TileSpacing', 'compact');
        for iT=1:nPlot
            nexttile
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,reshape(PlotUs(:,iT),Nx,Nx));
            title(strcat('$t=$',num2str(PlotTs(iT))))
            clim([min(PlotUs(:)) max(PlotUs(:))])
            %clim([0 3])
            set(gca,'YDir','Normal')
            colormap(turbo)
            hold on
            if (~isempty(Xf))
                plot(PlotXfs{iT}(:,1),PlotXfs{iT}(:,2),'ko','MarkerSize',0.2)
            end
            if (iT==nPlot)
                colorbar
            end
            if (iT==1)
                ylabel('$y$')
            end
            xlabel('$x$')
        end
        nexttile
        imagesc(rSim,tSim,XCorsSim/max(abs(XCorsSim(:))))
        clim([-1 1])
        colorbar
        xlabel('$\Delta r$')
        ylabel('$\Delta t$')
        title('Rho-actin Xcor')
    end
    if (0)
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
        imagesc(Uvals,dtvals,DistsByR/max(abs(DistsByR(:))))
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
        %close all;
    end
end
