function [Xf,nPerFil,TimeAtMax] = SetUpActin(L,ds,Xf,PoreSize,nFil,Lf)
    % Set up actin network
    if (isempty(Xf))
        if (isempty(PoreSize) || PoreSize<=0)
            Xf=[];
            for iFil=1:nFil
                x0=rand(1,2)*L;
                tau=randn(1,2);
                tau=tau/norm(tau);
                xf=x0+(0:ds:Lf)'.*tau;
                Xf=[Xf; xf];
            end
            nPerFil=length(0:ds:Lf)*ones(nFil,1);
            TimeAtMax = zeros(nFil,1);
        else
            % Regular grid with given pore size
            nFil = L/PoreSize;
            Lf=L;
            for iP=0:nFil-1
                x0 = [0 PoreSize*iP];
                tau=[1 0];
                xf=x0+(0:ds:Lf-ds)'.*tau;
                Xf=[Xf; xf];
                x0 = [PoreSize*iP 0];
                tau=[0 1];
                xf=x0+(0:ds:Lf-ds)'.*tau;
                Xf=[Xf; xf];
            end
            nFil=2*nFil;
            nPerFil=length(0:ds:Lf-ds)*ones(nFil,1);
            TimeAtMax = zeros(nFil,1);
        end
    end
end