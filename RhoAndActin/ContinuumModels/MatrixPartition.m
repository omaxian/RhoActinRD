function Regions = MatrixPartition(Nreg,L,Nx,dx)
    % Random pts;
    pts = rand(Nreg,2).*[L L];
    % All pts in the matrix that are closer to those points (respecting
    % periodicity)
    Regions = zeros(Nx);
    for iX=1:Nx
        for iY=1:Nx
            mpt = [(iX-1)*dx (iY-1)*dx];
            dists = mpt-pts;
            dists = dists - round(dists/L)*L;
            ndists = sqrt(sum(dists.*dists,2));
            [~,minpt]=min(ndists);
            Regions(iY,iX)=minpt;
        end
    end
end
