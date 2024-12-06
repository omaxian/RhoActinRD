function FiltData = FilterData(Data,nModes,nModesTime)  
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
    % Filter in time
    if (nModesTime > 0)
        for iX=1:nx
            for iY=1:ny
                TimeCoarse=reshape(FiltData(iX,iY,:),1,nFr);
                pTime = [TimeCoarse fliplr(TimeCoarse(2:end-1))];
                kvv=[0:nFr-1 -nFr+2:-1];
                pHat = fft(pTime);
                pHat(abs(kvv)>nModesTime)=0;
                FiltTime = ifft(pHat);
                FiltData(iX,iY,:)=FiltTime(1:nFr);
            end
        end
    end
end