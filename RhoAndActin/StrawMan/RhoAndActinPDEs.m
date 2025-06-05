function [Statistics,st] = RhoAndActinPDEs(Params,dt,postproc,RandomNuc,...
    seed)
    rng(seed);
    MakeMovie=1;
    koff0=Params(1);
    rf = Params(2);
    Nuc0=Params(3);
    NucEn=Params(4);
    koffAct = Params(5);
    Du = 0.1; % The size of the waves depends on Du
    Dv = Params(6);
    kbasal=Params(7);
    kfb=Params(8);
    KFB=Params(9);
    ChangeEvery = 1/koffAct;
    L=20;
    Nx=100; 
    % Solve for the steady states 
    [rts,rtstab,mmg,nmg] = PDERoots(Params,Du,L,Nx);
    if (isscalar(rts(:,1)))
        ICRange = [0 2*rts(:,1); 0*rts(:,2) 2*rts(:,2)];
    else 
        ICRange = [min(rts(:,1)) max(rts(:,1)); ...
            min(rts(:,2)) max(rts(:,2))];
    end
    
    %dt = 0.1; % Stability limit is 1
    tf = 220;
    tsaves = [40];
    saveEvery=floor(1e-6+1/dt);
    ChangeEverySteps = ceil(ChangeEvery/dt);
    
    % Initialize grid
    dx=L/Nx;
    x=(0:Nx-1)*dx;
    y=(0:Nx-1)*dx;
    [xg,yg]=meshgrid(x,y);
    % Set up a wave initial condition (works for worms)
    if (Dv>Du)
    u = ICRange(1,1)+0.5*(1+sin(mmg*2*pi*xg/L).*sin(nmg*2*pi*yg/L))...
       *(ICRange(1,2)-ICRange(1,1));
    v = ICRange(2,1)+0.5*(1+sin(mmg*2*pi*xg/L).*sin(nmg*2*pi*yg/L))...
       *(ICRange(2,2)-ICRange(2,1)); 
    else
    load('Initials.mat') %(starfish)
    u = uIC;
    v = sum(vIC,3)*2*pi/8;
    end
    kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
    [kx,ky]=meshgrid(kvals);
    ksq=kx.^2+ky.^2;
    DivFacFourier_u = (1/dt+Du*ksq);
    DivFacFourier_v = (1/dt+Dv*ksq);
    nSt = floor(tf/dt+1e-10);
    nSave = floor(tf/(saveEvery*dt)+1e-10);
    AllRho=zeros(Nx,Nx,nSave);
    AllActin=zeros(Nx,Nx,nSave);

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
        if (RandomNuc && mod(iT-1,ChangeEverySteps)==0)
            % Each spatial region gets a random theta for nucleation
            Nuc0s = zeros(Nx);
            NucEns = zeros(Nx);
            % % Make rates uniform in theta but not x
            SquareRegionSize = 4; % in um^2
            Nreg = ceil(L^2/SquareRegionSize);
            Regions = MatrixPartition(Nreg,L,Nx,dx);
            % For spatial anisotropy
            Vals0ByRegion = 2*rand(Nreg,1)*Nuc0;
            Vals1ByRegion = 2*rand(Nreg,1)*NucEn;
            for iX=1:Nx
                for iY=1:Nx
                    Nuc0s(iY,iX,:)=Vals0ByRegion(Regions(iY,iX));
                    NucEns(iY,iX,:)=Vals1ByRegion(Regions(iY,iX));
                end
            end
        elseif (~RandomNuc)
            Nuc0s=Nuc0*ones(Nx);
            NucEns=NucEn*ones(Nx);
        end
        RHS_u = (kbasal+kfb*u.^3./(KFB+u.^3))-(koff0+rf*v).*u;
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
        if (sum(abs(iT*dt-tsaves)<1e-10)>0)
            index = find(abs(iT*dt-tsaves)<1e-10);
            PlotUs(:,index)=reshape(u,[],1);
            PlotVs(:,index)=reshape(v,[],1);
            PlotTs(index)=iT*dt;
        end
        if (mod(iT,50)==0 && MakeMovie)
            tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
            nexttile
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
            set(gca,'YDir','Normal')
            %colorbar
            clim(ICRange(1,:))
            %colormap("turbo")
            title(sprintf('$t= %1.1f$',iT*dt))
            %clim(ICRange(1,:))
            colormap(turbo)
            colorbar
            pbaspect([1 1 1])
            nexttile
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,v);
            set(gca,'YDir','Normal')
            pbaspect([1 1 1])
            colormap(turbo)
            clim(ICRange(2,:))
            colorbar
            movieframes(iT)=getframe(f);
        end
        if (mod(iT-1,saveEvery)==0)
            AllRho(:,:,(iT-1)/saveEvery+1)=u;
            AllActin(:,:,(iT-1)/saveEvery+1)=v;
        end
    end
    % Post-process to get cross correlations and excitation sizes
    % Compute cross correlation function
    if (postproc)
        AllActin=AllActin(:,:,21:end);
        AllRho=AllRho(:,:,21:end);
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
    KymoPlot
    keyboard
end