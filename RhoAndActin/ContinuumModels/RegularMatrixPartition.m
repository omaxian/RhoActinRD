function Regions = RegularMatrixPartition(Nreg1D,L,Nx)
    % Random pts;
    dxR = L/Nreg1D;
    % All pts in the matrix that are closer to those points (respecting
    % periodicity)
    Regions = zeros(Nx);
    for iX=1:Nx
        for iY=1:Nx
            mpt = [(iX-1) (iY-1)]*L/Nx;
            regNum = floor(mpt/dxR);
            regNumAll = regNum(1)*Nreg1D + regNum(2) + 1;
            Regions(iY,iX)=regNumAll;
        end
    end
end
