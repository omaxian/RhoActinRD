function [rts,stability,minM,minN] = PDERoots(Params,Du,L,Nx)
    koff0=Params(1);
    rf = Params(2);
    Nuc0=Params(3);
    NucEn=Params(4);
    koffAct = Params(5);
    Dv=Params(6);
    kbasal=Params(7);
    kfb=Params(8);
    KFB=Params(9);
    % Coefficients if v_t = (n0+NucEn*u)-koffAct*v
    % pcoeffs=[-NucEn*rf/koffAct; -koff0-Nuc0*rf/koffAct;...
    %    kbasal+kfb; -KFB*NucEn*rf/koffAct; -KFB*koff0-KFB*Nuc0*rf/koffAct;...
    %    kbasal*KFB];
    % Coefficients if v_t = (n0+NucEn*u^2)-koffAct*v
    pcoeffs=[-NucEn*rf/koffAct; 0; -koff0-Nuc0*rf/koffAct;...
       kbasal+kfb-KFB*NucEn*rf/koffAct; 0; -KFB*koff0-KFB*Nuc0*rf/koffAct;...
       kbasal*KFB];
    urts = roots(pcoeffs);
    % Throw out complex and negative roots
    urts(imag(urts)~=0)=[];
    urts(urts<0)=[];
    vrts= (Nuc0+NucEn*urts.^2)/koffAct;
    urts(vrts<0)=[];
    vrts(vrts<0)=[];
    rts=[urts vrts];
    rts=unique(rts,'rows');
    Nrt=length(urts);
    stability=ones(Nrt,1);
    J11 = 3*kfb*KFB*urts.^2./(KFB+urts.^3).^2 - (koff0+rf*vrts);
    J12 = -rf*urts;
    J21 = 2*NucEn*urts;
    J22 = -koffAct*ones(length(urts),1);
    minM=0;
    minN=0;
    maxEig=-inf;
    detJs = J11.*J22-J21.*J12;
    trJs = J11+J22;
    %Stable = detJs > 0 & trJs < 0;
    Unstable = detJs < 0 | (detJs > 0 & trJs > 0);
    %spiral = 4*detJs - trJs.^2 > 0;
    stability(Unstable)=-1;
    if (Du > 0 || Dv > 0)
        for m=1:Nx-1
            for n=1:Nx-1
                kxsq = 4*pi^2*m^2/L^2;
                kysq = 4*pi^2*n^2/L^2;
                J11D = J11 - (kxsq+kysq)*Du;
                J22D = J22 - (kxsq+kysq)*Dv;
                for iP=1:length(J11D)
                maxEigMN = max(real(eig([J11D(iP) J12(iP); J21(iP) J22D(iP)])));
                if (m+n>0 && maxEigMN > maxEig)
                    maxEig=maxEigMN;
                    minM=m;
                    minN=n;
                end
                end
            end
        end
    end
end

