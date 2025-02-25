function [rts,stability,spiral] = findFNRoots(a,b,Iext,tau,Du,Dv,Nx)
    if (b==0)
        urts=-a;
    else
        pcoeffs = [1/3;0;1/b-1;a/b-Iext];
        urts = roots(pcoeffs);
        % Throw out complex roots
        urts(imag(urts)~=0)=[];
    end
%     us=-5:0.1:5;
%     v1=us-1/3*us.^3+Iext;
%     v2=1/b*(us+a);
%     plot(us,v1)
%     hold on
%     plot(us,v2)
    vrts= urts - 1/3*urts.^3+Iext;
    rts=[urts vrts];
    % Get stability of the roots
    Nrt=length(urts);
    stability=zeros(Nrt,1);
    detJs = (b*urts.^2-b+1)/tau;
    trJs = 1-urts.^2 - b/tau;
    Stable = detJs > 0 & trJs < 0;
    Unstable = detJs < 0 | (detJs > 0 & trJs > 0);
    spiral = 4*detJs - trJs.^2 > 0;
    if (Du > 0 || Dv > 0)
        for m=1:Nx
            ksq = 8*pi^2*m^2;
            detJs = (1-urts.^2-ksq*Du)*(-b/tau-ksq*Dv)+1/tau;
            trJs = (1-urts.^2-ksq*Du)+(-b/tau-ksq*Dv);
            Stablem = detJs > 0 & trJs < 0;
            Unstablem = detJs < 0 | (detJs > 0 & trJs > 0);
            spiralm = 4*detJs - trJs.^2 > 0;
            if (Unstablem)
                sprintf('Mode %d unstable',m)
            end
        end
    end
    stability(Stable)=1;
    stability(Unstable)=-1;
end

