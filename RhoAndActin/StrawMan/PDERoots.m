function [rts,stability,spiral] = PDERoots(Params)
    kbasal=Params(6);
    kfb=Params(7);
    KFB=Params(8);
    %Du=0.1;
    koff0=Params(1);
    rf = Params(2);
    Nuc0=Params(3);
    NucEn=Params(4);
    %Dv=Params(5);
    koffAct = Params(5);
    pcoeffs=[-NucEn*rf/koffAct; -koff0-Nuc0*rf/koffAct;...
       kbasal+kfb; -KFB*NucEn*rf/koffAct; -KFB*koff0-KFB*Nuc0*rf/koffAct;...
       kbasal*KFB];
    urts = roots(pcoeffs);
    % Throw out complex and negative roots
    urts(imag(urts)~=0)=[];
    urts(urts<0)=[];
    vrts= (Nuc0+NucEn*urts)/koffAct;
    urts(vrts<0)=[];
    vrts(vrts<0)=[];
    rts=[urts vrts];
    % u1=max(min(urts)-0.5,0):0.05:max(urts)+2;
    % v1=max(min(vrts)-0.5,0):0.02:max(vrts)+0.5;
    % figure(3)
    % [u,v]=meshgrid(u1,v1);
    % RHS_u = (kbasal+kfb*u.^3./(KFB+u.^3))-(koff0+rf*v).*u;
    % RHS_v = (Nuc0+NucEn*u.*v) - koffAct*v;
    % plot(u1,1/rf*(-koff0+1./u1.*(kbasal+kfb*u1.^3./(KFB+u1.^3))))
    % hold on
    % plot(u1,1/koffAct*(Nuc0+NucEn*u1))
    % quiver(u,v,RHS_u,RHS_v)
    % hold on
    % Get stability of the roots
    Nrt=length(urts);
    stability=zeros(Nrt,1);
    J11 = 3*kfb*KFB*urts.^2./(KFB+urts.^3).^2 - (koff0+rf*vrts);
    J12 = -rf*urts;
    J21 = NucEn;
    J22 = -koffAct;
    detJs = J11.*J22-J21.*J12;
    trJs = J11+J22;
    Stable = detJs > 0 & trJs < 0;
    Unstable = detJs < 0 | (detJs > 0 & trJs > 0);
    spiral = 4*detJs - trJs.^2 > 0;
    % if (Du > 0 || Dv > 0)
    %     for m=1:Nx
    %         ksq = 8*pi^2*m^2;
    %         detJs = (1-urts.^2-ksq*Du)*(-b/tau-ksq*Dv)+1/tau;
    %         trJs = (1-urts.^2-ksq*Du)+(-b/tau-ksq*Dv);
    %         Stablem = detJs > 0 & trJs < 0;
    %         Unstablem = detJs < 0 | (detJs > 0 & trJs > 0);
    %         spiralm = 4*detJs - trJs.^2 > 0;
    %         if (Unstablem)
    %             sprintf('Mode %d unstable',m)
    %         end
    %     end
    % end
    stability(Stable)=1;
    stability(Unstable)=-1;
end

