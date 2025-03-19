function [Uvals,dtvals,DistsByR] = CrossCorrelations(dx,dy,dt,Adata,Bdata,padxy)
    Adata=Adata-mean(Adata(:));
    Bdata=Bdata-mean(Bdata(:));
    [Ny,Nx,Nt]=size(Adata);
    dtvals=(-Nt+1:Nt-1)*dt;
    if (padxy)
        az = zeros(2*Ny-1,2*Nx-1,2*Nt-1);
        bz = zeros(2*Ny-1,2*Nx-1,2*Nt-1);
        dxvals=(-Nx+1:Nx-1)*dx;
        dyvals=(-Ny+1:Ny-1)*dy;
        shft = [Ny-1 Nx-1 Nt-1];
    else
        az = zeros(Ny,Nx,2*Nt-1);
        bz = zeros(Ny,Nx,2*Nt-1);
        dxvals=(-Nx/2:Nx/2-1)*dx;
        dyvals=(-Ny/2:Ny/2-1)*dy;
        shft = [Ny/2 Nx/2 Nt-1];
    end
    az(1:Ny,1:Nx,1:Nt)=Adata;
    bz(1:Ny,1:Nx,1:Nt)=Bdata;
    c3=circshift(ifftn(fftn(bz).*conj(fftn(az))),shft);
    c3_2D=reshape(c3,size(c3,1)*size(c3,2),size(c3,3));
    
    [xg,yg]=meshgrid(dxvals,dyvals);
    rg=sqrt(xg(:,:,1).^2+yg(:,:,1).^2);
    [val,SortedInds]=sort(rg(:));
    [Uvals,~,ic]=uniquetol(val,1e-2);
    % The map is SortedInds -> ic
    nu=length(Uvals);
    DistsByR = zeros(length(dtvals),nu);
    if (Nx==100 & dx==0.2)
        load 'CrossCorrMapFromXYToR.mat' Inds
        for iP=1:nu
            DistsByR(:,iP)=mean(c3_2D(Inds{iP},:),1);
        end
    else
    warning('Should not be here!')
    for iT=1:length(dtvals)
        ThisImg = c3(:,:,iT);
        for iP=1:nu
            DistsByR(iT,iP)=mean(ThisImg(SortedInds(ic==iP)));
        end
    end 
    end
end



