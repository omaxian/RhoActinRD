function [Uvals,dtvals,DistsByR] = CrossCorrelations(dx,dy,dt,Adata,Bdata)
    Adata=Adata-mean(Adata(:));
    Bdata=Bdata-mean(Bdata(:));
    [Ny,Nx,Nt]=size(Adata);
    % Pad in the z dimension
    az = Adata;
    az(:,:,end:end+Nt-1)=0;
    bz = Bdata;
    bz(:,:,end:end+Nt-1)=0;
    c3=circshift(ifftn(fftn(bz).*conj(fftn(az))),[Ny/2 Nx/2 Nt-1]);
    dtvals=[(-Nt+1:-1)*dt (0:Nt-1)*dt];
    dxvals=(-Nx/2:Nx/2-1)*dx;
    dyvals=(-Ny/2:Ny/2-1)*dy;
    [xg,yg]=meshgrid(dxvals,dyvals);
    rg=sqrt(xg(:,:,1).^2+yg(:,:,1).^2);
    [val,SortedInds]=sort(rg(:));
    [Uvals,~,ic]=uniquetol(val,1e-2);
    % The map is SortedInds -> ic
    nu=length(Uvals);
    DistsByR = zeros(length(dtvals),nu);
    for iT=1:length(dtvals)
        ThisImg = c3(:,:,iT);
        for iP=1:nu
            DistsByR(iT,iP)=mean(ThisImg(SortedInds(ic==iP)));
        end
    end
    figure
    title('Cross correlation')
    imagesc(Uvals,dtvals,DistsByR)
    xlabel('$\Delta r$')
    ylabel('$\Delta t$')
end



