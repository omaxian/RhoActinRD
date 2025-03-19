function [Stats,unstable] = FNDynamics(Params,seed,dt,MakeMovie)
    rng(seed);
    L = 20;
    a = Params(1);
    b = Params(2);
    Iext = Params(3);
    tau = Params(4);
    Du = 0.1; % The size of the waves depends on Du
    Dv = Params(5);
    tf = 250;
    unstable = 0;
    
    % Numerical params
    %dt = 0.5;
    saveEvery=floor(1e-6+1/dt);
    Nx = 400; % The grid spacing has GOT to resolve the boundary layer - scales with min(D)!
    dx = L/Nx;
    kvals = [0:Nx/2 -Nx/2+1:-1]*2*pi/L;
    [kx,ky]=meshgrid(kvals);
    ksq=kx.^2+ky.^2;
    DivFacFourier_u = (1/dt+Du*ksq);
    DivFacFourier_v = (1/dt+Dv*ksq);
    nSave = floor(tf/(saveEvery*dt)+1e-10);
    AllRho=zeros(Nx,Nx,nSave);
    AllActin=zeros(Nx,Nx,nSave);
    
    % Steady states
    [rts,~,~] = findFNRoots(a,b,Iext,tau,Du,Dv,Nx);
    [nRts,~]=size(rts);
    rts=real(rts);
    if (nRts>1)
        RtRange = max(rts)-min(rts);
        RtStart = min(rts)-1.1*RtRange;
        RtEnd = max(rts)+1.1*RtRange;
    else
        RtRange = 1;
        RtStart = rts-0.5*RtRange;
        RtEnd = rts+0.5*RtRange;
    end
    u = RtStart(1)+(RtEnd(1)-RtStart(1)).*rand(Nx);
    v = RtStart(2)+(RtEnd(2)-RtStart(2)).*rand(Nx);
    dx=L/Nx;
    x=(0:Nx-1)*dx;
    y=(0:Nx-1)*dx;
    [xg,yg]=meshgrid(x,y);
    u = RtStart(1)+0.5*(1+sin(2*pi*xg/L).*...
    sin(2*pi*yg/L))*(RtEnd(1)-RtStart(1));
    v = RtStart(2)+0.5*(1+sin(2*pi*xg/L).*...
    sin(2*pi*yg/L))*(RtEnd(2)-RtStart(2));
    
    if (MakeMovie)
        f=figure('Position',[100 100 600 300]);
        %tiledlayout(1,4,'Padding', 'none', 'TileSpacing', 'compact');
    end
    nSt = tf/dt;
    for iT=1:nSt
        RHS_u = u-1/3*u.^3-v+Iext;
        RHS_v = 1/tau*(u+a-b*v);
        RHSHat_u = fft2(u/dt+RHS_u);
        uHatNew = RHSHat_u./(DivFacFourier_u);
        u = ifft2(uHatNew);
        RHSHat_v = fft2(v/dt+RHS_v);
        vHatNew = RHSHat_v./(DivFacFourier_v);
        v = ifft2(vHatNew);
        % Print error if it goes unstable
        if (max(abs(u(:)))>1e5)
            %warning('Unstable simulation - returning')
            Stats.XCor=0;
            Stats.MeanActin=inf;
            unstable=1;
            return
        end
        if (mod(iT-1,10)==0 && MakeMovie)
            tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
            nexttile
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
            colormap turbo
            rngU=range(u(:));
            clim([min(u(:)) max(u(:))])
            %clim([-2.1 3.6])
            if (rngU<1)
                clim([mean(u(:))-0.5 mean(u(:))+0.5])
            end
            xlabel('$x$ ($\mu$m)')
            ylabel('$y$ ($\mu$m)')
            %colorbar
            pbaspect([1 1 1])
            title(strcat('$\rho$, $t=$',num2str(iT*dt)))
            nexttile
            imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,v);
            colormap turbo
            rngV=range(v(:));
            clim([min(v(:)) max(v(:))])
            %clim([1.2 2.2])
            xlabel('$x$ ($\mu$m)')
            if (rngV<1)
                clim([mean(v(:))-0.5 mean(v(:))+0.5])
            end
            %colorbar
            title(strcat('$f$, $t=$',num2str(iT*dt)))
            pbaspect([1 1 1])
            movieframes((iT-1)/10+1)=getframe(f);
        end
        % if (mod(iT,30)==0 && MakeMovie && iT*dt>100)
        %     nexttile
        %     imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,u);
        %     colormap turbo
        %     rngU=range(u(:));
        %     clim([min(u(:)) max(u(:))])
        %     if (rngU<1)
        %         clim([mean(u(:))-0.5 mean(u(:))+0.5])
        %     end
        %     colorbar
        %     clim([-2 2])
        %     pbaspect([1 1 1])
        %     title(strcat('$u, t=$',num2str(iT*dt)))
        % end
        if (mod(iT-1,saveEvery)==0)
            AllRho(:,:,(iT-1)/saveEvery+1) = u;
            AllActin(:,:,(iT-1)/saveEvery+1) = v;
            % Return if the range in space and time is zero
            if (iT>1)
                RhosBeforeAndNow=[AllRho(:,:,(iT-1)/saveEvery+1); ...
                    AllRho(:,:,(iT-1)/saveEvery,:)];
                if (range(RhosBeforeAndNow)<1e-4)
                    % Break out 
                    Stats.XCor=0;
                    ActinSoFar=AllActin(:,:,1:(iT-1)/saveEvery+1);
                    Stats.MeanActin=mean(ActinSoFar(:));
                    return;
                end
            end
        end
    end
    % Post-process to get cross correlations and excitation sizes
    % Compute cross correlation function
    AllActin=AllActin(:,:,51:end);
    AllRho=AllRho(:,:,51:end);
    Thres=mean(AllRho(:));
    if (length(rts(:,1))>1)
        % Take the middle root
        SortedRoots=sort(rts(:,1),'ascend');
        Thres1=SortedRoots(2);
    else
        Thres1=Thres;
    end
    RhoThres=AllRho>Thres;
    RhoThres1=AllRho>Thres1;
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
    ExSizes1=[];
    for iT=1:nFr
        CC = bwconncomp(RhoThres1(:,:,iT));
        L2=CC2periodic(CC,[1 1],'L');
        NumExcitations(iT)=max(L2(:));
        for iJ=1:max(L2(:))
            ExSizes1=[ExSizes1;sum(L2(:)==iJ)/(Nx^2)*L^2];
        end
    end
    [rSim,tSim,XCorsSim] = CrossCorrelations(dx,dx,dt*saveEvery,...
        AllRho,AllActin,0);
    Stats.XCor=XCorsSim/max(abs(XCorsSim(:)));
    Stats.rSim=rSim;
    Stats.tSim=tSim;
    Stats.ExSizes=ExSizes;
    Stats.ExSizesStThres=ExSizes1;
    Stats.NumExcitations=NumExcitations;
    Stats.MeanActin=mean(AllActin(:));
end

% Bds=[floor(10*min(uv1)-3); ceil(10*max(uv1)+3)]/10;
% %Bds=[-3 -3; 3 3];
% us=Bds(1,1):0.1:Bds(2,1);
% vs=Bds(1,2):0.1:Bds(2,2);
% [us,vs]=meshgrid(us,vs);
% dudt = us-1/3*us.^3-vs+Iext;
% dvdt = 1/tau*(us+a-b*vs);
% quiver(us,vs,dudt,dvdt)
% hold on
% plot(uv1(:,1),uv1(:,2))
% xlim(Bds(:,1)')
% ylim(Bds(:,2)')
% plot(rts(:,1),rts(:,2),'ko','MarkerFaceColor','k')
% %imagesc((0:Nx-1)*dx,(0:Nx-1)*dx,reshape(uv(:,2),Nx,Nx));
% %clim([-0.3411    1.2060])
% %title(strcat('$D_u=$',num2str(Du),' $D_v=$',num2str(Dv)))