function [Statistics,st] = RhoAndActinPDEs(Params,dt,seed,doPlot)
    %rng(seed);
    MakeMovie=0;
    kbasal=Params(7);
    kfb=Params(8);
    KFB=Params(9);
    koff0=Params(1);
    rf = Params(2);
    Nuc0=Params(3);
    NucEn=Params(4);
    koffAct = Params(5);
    % Solve for the steady states 
    [rts,~,~] = PDERoots([Params(1:5);Params(7:9)]);
    if (isscalar(rts(:,1)))
        ICRange = [0 2*rts(:,1); 0 2*rts(:,2)];
    elseif (length(rts(:,1))==3)
        ICRange = [min(rts(:,1)) max(rts(:,1)); ...
            min(rts(:,2)) max(rts(:,2))];
    end
    
    %dt = 0.1; % Stability limit is 1
    tf = 200;
    Du = 0.1; % The size of the waves depends on Du
    Dv = Params(6);
    tsaves = [160];
    saveEvery=floor(1e-6+1/dt);
    
    % Initialize grid
    L=20;
    Nx=100; % The grid spacing 
    dx=L/Nx;
    x=(0:Nx-1)*dx;
    y=(0:Nx-1)*dx;
    [xg,yg]=meshgrid(x,y);
    % Set up a wave initial condition
    %u = ICRange(1,1)+0.5*(1+sin(2*pi*xg/L).*sin(2*pi*yg/L))*(ICRange(1,2)-ICRange(1,1));
    %v = ICRange(2,1)+0.5*(1+sin(2*pi*xg/L).*sin(2*pi*yg/L))*(ICRange(2,2)-ICRange(2,1));
    u = ICRange(1,1)+rand(Nx)*(ICRange(1,2)-ICRange(1,1));
    v = ICRange(2,1)+rand(Nx)*(ICRange(2,2)-ICRange(2,1)); 
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
        % Break out of sims not doing anything
        if (range(u(:)) <1e-2 && range(v(:)) < 1e-2)
            Statistics.XCor=0;
	        Statistics.MeanActin=0;
            return;
        end
        RHS_u = (kbasal+kfb*u.^3./(KFB+u.^3))-(koff0+rf*v).*u;
        RHS_v = (Nuc0+NucEn*u) - koffAct*v;
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
        if (mod(iT,20)==0 && MakeMovie)
            tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
            nexttile
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
            set(gca,'YDir','Normal')
            %colorbar
            %clim([0 StSt(end)])
            %colormap("turbo")
            title(sprintf('$t= %1.1f$',iT*dt))
            clim(ICRange(1,:))
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
    AllActin=AllActin(:,:,41:end);
    AllRho=AllRho(:,:,41:end);
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
end