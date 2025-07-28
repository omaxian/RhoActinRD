% This is tma
function Statistics = RhoAndActinBasalNuc(Params,seed,doPlot)
    % Parameters:
    % koff0, rf, FullLifetime, Grow, Shrink, MaxLength, Nuc0, NucEn 
    % in that order
    % Output is the difference in the cross correlations compared to
    % experimental data
    rng(seed);
    MakeMovie=doPlot;
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
    Nuc0=Params(7);
    NucEn=Params(8)/max(StSt)^2;
    dt = 0.25; % Stability limit is 1
    tf = 241;
    Du=0.1; % The size of the waves depends on Du
    tsaves = [40 240];
    % Parameters for the actin
    PoreSize=[];
    ds=0.1;
    FullLifetime = Params(3); % the average lifetime of a filament (s)
    GrowRate=Params(4);
    ShrinkRate=Params(5);
    GrowAmt = (GrowRate/ds*dt);% (#mon per time step - first number is um/s)
    ShrinkAmt = (ShrinkRate/ds*dt); % (#mon per time step - first number is um/s)
    MaxLength = max(ds,floor(Params(6)/ds)*ds);
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
    [xg,yg]=meshgrid(x,y);
    u = ones(Nx,Nx)*ICScale;
    %u(xg/L<0.4 | xg/L> 0.6 | yg/L < 0.4 | yg/L > 0.6)=min(StSt);
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
        f=figure('Position',[100 100 500 500]);
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
            NucRate=(Nuc0+NucEn*u.^2)*dt; % # nucleates/length^2
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
        if (mod(iT-1,40)==0 && MakeMovie)
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
            set(gca,'YDir','Normal')
            hold on
            if (~isempty(Xf))
                Xfpl=Xf-floor(Xf/L)*L;
                scatter(Xfpl(:,1),Xfpl(:,2),10,'o','filled',...
                    'MarkerFaceColor',[0.87    0.49    0.0],...
                    'MarkerEdgeColor','None',...
                    'MarkerFaceAlpha',0.5)
            end
            %colorbar
            %clim([0 StSt(end)])
            %colormap("turbo")
            title(sprintf('$t= %1.1f$',t-40))
            clim([0 max(StSt)])
            colormap(sky)
            hold off
            pbaspect([1 1 1])
            xlabel('$x$ ($\mu$m)')
            ylabel('$y$ ($\mu$m)')
            yticks(0:5:15)
            xticks(0:5:15)
            movieframes(iT)=getframe(f);
        end
        if (mod(iT-1,saveEvery)==0)
            AllRho(:,:,(iT-1)/saveEvery+1)=u;
            AllActin(:,:,(iT-1)/saveEvery+1)=fg;
            % Get x coordinates of all filaments in middle
            Xfpl=Xf-floor(Xf/L)*L;
            kymopt=Nx/2;
            xCoords{(iT-1)/saveEvery+1}=...
                Xfpl(Xfpl(:,2)>(kymopt-1)*dx & Xfpl(:,2)<kymopt*dx,1);
        end
    end
    % Post-process to get cross correlations and excitation sizes
    % Compute cross correlation function
    BurnIn=40;
    AllActin=AllActin(:,:,BurnIn+1:end);
    AllRho=AllRho(:,:,BurnIn+1:end);
    xCoords=xCoords(BurnIn+1:end);
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
    Statistics.XCor=XCorsSim/max(abs(XCorsSim(:)));
    Statistics.rSim=rSim;
    Statistics.tSim=tSim;
    Statistics.ExSizes=ExSizes;
    Statistics.NumExcitations=NumExcitations;
    Statistics.MeanActin=mean(AllActin(:));
    if (0)
        % Snapshots
        [~,nPlot]=size(PlotUs);
        tiledlayout(1,nPlot,'Padding', 'none', 'TileSpacing', 'compact');
        for iT=1:2
            figure(iT+1)
            nexttile
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,reshape(PlotUs(:,iT),Nx,Nx));
            %title(strcat('$t=$',num2str(PlotTs(iT))))
            clim([min(AllRho(:)) max(AllRho(:))])
            %clim([0 3])
            set(gca,'YDir','Normal')
            colormap(sky)
            hold on
            pbaspect([1 1 1])
            if (~isempty(Xf))
                scatter(PlotXfs{iT}(:,1),PlotXfs{iT}(:,2),2.5,'o','filled',...
                    'MarkerFaceColor',[0.87    0.49    0],...
                    'MarkerEdgeColor','None',...
                    'MarkerFaceAlpha',0.35)
            end
            % if (iT==nPlot)
            %     colorbar
            % end
            % if (iT==1)
            %     ylabel('$y$ ($\mu$m)')
            % end
            %xlabel('$x$ ($\mu$m)')
            xticklabels('')
            yticklabels('')
        end
        % Kymograph
        figure(4)
        nexttile
        RhoT=reshape(AllRho(kymopt,:,:),Nx,tf-BurnIn)';
        imagesc((0:Nx-1)*dx,0:tf-BurnIn-1,RhoT)
        colormap sky
        hold on
        for j=1:tf-BurnIn
            xpl=xCoords{j};
            scatter(xpl,(j-1)*ones(length(xpl),1),2,'s', ...
                    'MarkerFaceColor',[0.87    0.49    0],...
                    'MarkerEdgeColor','None',...
                    'MarkerFaceAlpha',0.04)
        end
        %xlabel('$x$ ($\mu$m)')
        xticklabels('')
        pbaspect([1 1.25 1])
        yticklabels('')
    end
end