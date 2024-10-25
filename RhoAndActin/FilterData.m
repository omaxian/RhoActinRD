function FiltData = FilterData(Data,nModes)  
    [ny,nx,nFr]=size(Data);
    Nx=2*nx-2;
    FiltData=zeros(ny,nx,nFr);
    kvalsX = [0:nx-1 -nx+2:-1]*2*pi;
    kvalsY = [0:ny-1 -ny+2:-1]*2*pi;
    kFilt = (nModes*2*pi)^2;
    [kx,ky]=meshgrid(kvalsX,kvalsY);
    ksq=kx.^2+ky.^2;
    for iT=1:nFr
        ThisFr = Data(:,:,iT);
        PerExtData = [ThisFr fliplr(ThisFr(:,2:end-1))];
        PerExtData = [PerExtData; flipud(PerExtData(2:end-1,:))];
        DataHat = fft2(PerExtData);
        % Add zero wave numbers or subtract wave numbers
        % Zero unpaired mode
        DataHat(:,nx)=0;
        DataHat(ny,:)=0;
        if (Nx > 2*nx-2)
            DataHat = [DataHat(:,1:nx) zeros(2*ny-2,Nx-(2*nx-2)) DataHat(:,nx+1:end)];
        else
            DataHat = [DataHat(:,1:Nx/2+1) DataHat(:,end-Nx/2+2:end)];
            % Zero out unpaired mode
            DataHat(:,Nx/2+1)=0;
        end
        if (Nx > 2*ny-2)
            DataHat = [DataHat(1:ny,:); zeros(Nx-(2*ny-2),Nx); DataHat(ny+1:end,:)];
        else
            DataHat = [DataHat(1:Nx/2+1,:); DataHat(end-Nx/2+2:end,:)];
            DataHat(Nx/2+1,:)=0;
        end
        DataHat(ksq>kFilt)=0;
        PaddedRev = ifft2(DataHat);
        FiltData(:,:,iT) = PaddedRev(1:ny,1:nx,:);
    end
end