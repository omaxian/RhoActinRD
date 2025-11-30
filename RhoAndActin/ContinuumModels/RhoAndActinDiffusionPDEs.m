function [Statistics,st] = RhoAndActinDiffusionPDEs(Params,dt,postproc,RandomNuc,...
    seed,pltkymo)
    rng(seed);
    MakeMovie=0;
    koff0=Params(1);
    rf = Params(2);
    Nuc0=Params(3);
    NucEn=Params(4);
    koffAct = Params(5);
    Du = 0.1; % The size of the waves depends on Du
    SquareRegionSize = Params(10); % in um^2
    Dv = Params(6);
    kbasal=Params(7);
    kfb=Params(8);
    KFB=Params(9);
    ChangeEvery = 2/koffAct;
    L=20;
    Nx=100; 
    % Solve for the steady states 
    [rts,~,~,~] = PDERoots(Params(1:9),Du,L,Nx);
    StSt=rts(:,1);
    if (isscalar(StSt))
        Thres = 0.5*StSt(1);
    else 
        Thres = 0.5*(StSt(2,1)+StSt(3,1));
    end
    kThres = 0.5;
    % Set the new k so that there is a single st st
    kLo=0.5;
    pChk=Params;
    StStLo=StSt;
    while (length(StStLo)>1)
        kLo=kLo-0.1;
        pChk(1)=kLo;
        [rtsChk,~,~,~] = PDERoots(pChk,Du,L,Nx);
        StStLo=rtsChk(:,1);
    end
    
    %dt = 0.1; % Stability limit is 1
    tf = 540;
    saveEvery=floor(1e-6+1/dt);
    ChangeEverySteps = ceil(ChangeEvery/dt);
    Nreg1D = round(L/sqrt(SquareRegionSize));
    Regions = RegularMatrixPartition(Nreg1D,L,Nx);
    
    % Initialize grid
    dx=L/Nx;
    % Set up IC
    u = rts(1,1)*ones(Nx);
    v = rts(1,2)*ones(Nx);
    koffz=koff0;
    
    kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
    [kx,ky]=meshgrid(kvals);
    ksq=kx.^2+ky.^2;
    DivFacFourier_u = (1/dt+Du*ksq);
    DivFacFourier_v = (1/dt+Dv*ksq);
    nSt = floor(tf/dt+1e-10);
    nSave = floor(tf/(saveEvery*dt)+1e-10);
    AllRho=zeros(Nx,Nx,nSave);
    AllActin=zeros(Nx,Nx,nSave);
    AllRhoHat=zeros(Nx,Nx,nSave);
    AllActinHat=zeros(Nx,Nx,nSave);
    Stimulated=zeros(nSave,1);

    if (MakeMovie)
        close all;
        f=figure('Position',[100 100 700 400]);
    end
    
    % Run simulation
    st=1;
    for iT=1:nSt
        if (RandomNuc && mod(iT-1,ChangeEverySteps)==0)
            % Each spatial region gets a random theta for nucleation
            Nuc0s = zeros(Nx);
            NucEns = zeros(Nx);
            % For spatial anisotropy
            ValsByRegion = 2*rand(Nreg1D^2,1);
            for iX=1:Nx
                for iY=1:Nx
                    Nuc0s(iY,iX)=ValsByRegion(Regions(iY,iX))*Nuc0;
                    NucEns(iY,iX)=ValsByRegion(Regions(iY,iX))*NucEn;
                end
            end
        elseif (~RandomNuc)
            Nuc0s=Nuc0*ones(Nx);
            NucEns=NucEn*ones(Nx);
        end
        if (length(StSt)>1)
            if (max(u(:)) > Thres)
                koffz=koff0;
            elseif (min(koffz(:))>kThres)
                koffz = kLo;
            end
        end
        RHS_u = (kbasal+kfb*u.^3./(KFB+u.^3))-(koffz+rf*v).*u;
        RHS_v = (Nuc0s+NucEns.*u.^2) - koffAct*v;
        RHSHat_u = fft2(u/dt+RHS_u);
        uHatNew = RHSHat_u./DivFacFourier_u;
        uNew = ifft2(uHatNew);
        RHSHat_v = fft2(v/dt+RHS_v);
        vHatNew = RHSHat_v./DivFacFourier_v;
        vNew = ifft2(vHatNew);
        u = uNew;
        v = vNew;
        if (max(abs(u(:)) > 1e5))
            st=0;
            warning('Rejecting because of unstable simulation')
            Statistics.XCor=0;
	        Statistics.MeanActin=0;
            return;
        end
        if (mod(iT,200)==0 && MakeMovie)
            tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
            ax1=nexttile;
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
            set(gca,'YDir','Normal')
            clim([0 max(rts(:,1))])
            colormap(ax1,sky);
            title(sprintf('Rho; $t= %1.1f$',iT*dt))
            pbaspect([1 1 1])
            xlabel("$x$ ($\mu$m)")
            ylabel("$y$ ($\mu$m)")
            colorbar
            ax2=nexttile;
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,v);
            set(gca,'YDir','Normal')
            pbaspect([1 1 1])
            C2=[0.87 0.49 0];
            C1=[0.95 0.9 0.9];
            Cmap=C1+(0:100)'/100.*(C2-C1);
            colormap(ax2,Cmap)
            clim([0 max(rts(:,2))])
            colorbar
            title(sprintf('F-actin; $t= %1.1f$',iT*dt))
            xlabel("$x$ ($\mu$m)")
            movieframes(iT)=getframe(f);
        end
        if (mod(iT-1,saveEvery)==0)
            AllRho(:,:,(iT-1)/saveEvery+1)=u;
            AllActin(:,:,(iT-1)/saveEvery+1)=v;
            AllRhoHat(:,:,(iT-1)/saveEvery+1)=fft2(u);
            AllActinHat(:,:,(iT-1)/saveEvery+1)=fft2(v);
            Stimulated((iT-1)/saveEvery+1)=min(koffz(:))<kThres;
        end
    end
    % Post-process to get cross correlations and excitation sizes
    % Compute cross correlation function
    if (postproc)
        if (length(StSt)>1)
            % Find the longest continuous interval where there was no
            % stimulation
            CCStim=bwconncomp(~Stimulated);
            p = regionprops(CCStim, 'Area'); 
            longest = 0;
            if (~isempty(p))
                [longest,ind] = max([p.Area]); % Find the maximum area
                IncludeMe = CCStim.PixelIdxList{ind};
            end
            if (longest<120/(dt*saveEvery))
                IncludeMe=40/(dt*saveEvery)+1:size(AllActin,3);
                EnoughExcitation=0;
            else
                EnoughExcitation=1;
            end
        else
            longest = 0;
            EnoughExcitation=1;
            IncludeMe=40/(dt*saveEvery)+1:size(AllActin,3);
        end
        % Remove initial transients
        AllActin=AllActin(:,:,IncludeMe);
        AllRho=AllRho(:,:,IncludeMe);
        nFour=10;
        AllActinHat=AllActinHat(1:nFour,1:nFour,IncludeMe);
        AllRhoHat=AllRhoHat(1:nFour,1:nFour,IncludeMe);
        % Save magnitude and autocorrelation of Fourier modes
        MeanActinHat = mean(abs(AllActinHat),3);
        MeanRhoHat = mean(abs(AllRhoHat),3);
        % Compute autocorrelations at times 0.5, 2, 5, and 10
        TimeAcor=max([0.5 2 4 6 12],saveEvery*dt);
        Lags = TimeAcor/(saveEvery*dt);
        ACorsRho = zeros(nFour,nFour,length(Lags));
        ACorsAct = zeros(nFour,nFour,length(Lags));
        for j=1:nFour
            for k=1:nFour
                tser = reshape(abs(AllRhoHat(j,k,:)),[],1);
                aCors = autocorr(tser,NumLags=max(Lags));
                tserA = reshape(abs(AllActinHat(j,k,:)),[],1);
                aCorsA = autocorr(tserA,NumLags=max(Lags));
                for iL=1:length(Lags)
                    ACorsRho(j,k,:)=aCors(Lags+1);
                    ACorsAct(j,k,:)=aCorsA(Lags+1);
                end
            end
        end
        
        % Check if the Rho concentration is always below/above the saddle pt
        nT = size(AllRho,3);
        Thres=StSt(2);
        AboveByTime = reshape(sum(sum(AllRho>Thres,1),2),nT,1);
        AboveByTime = AboveByTime/Nx^2;
        [rSim,tSim,XCorsSim] = CrossCorrelations(dx,dx,dt*saveEvery,...
            AllRho,AllActin,0);
        XCorsSim=XCorsSim/max(abs(XCorsSim(:)));
        ResampledT = -120:2:120;
        ResampledX = 0:0.5:10;
        InterpolatedSim=ResampleXCor(XCorsSim,tSim,rSim,...
                 ResampledX,ResampledT,11,121);
        Statistics.XCor=InterpolatedSim;
        Statistics.rSim=ResampledX;
        Statistics.tSim=ResampledT;
        Statistics.PctAboveSaddle = AboveByTime;
        Statistics.MeanRhoHat = MeanRhoHat(:);
        Statistics.MeanActinHat = MeanActinHat(:);
        Statistics.ACorsRho = ACorsRho(:);
        Statistics.ACorsAct = ACorsAct(:);
        Statistics.TimeACor = TimeAcor;
        Statistics.EnoughExcitation=EnoughExcitation;
        Statistics.LongestNoStim = longest*dt*saveEvery;
    else
        Statistics=[];
    end
    if (pltkymo)
        if (MakeMovie)
            close(f)
        end
        KymoPlot
    end
end