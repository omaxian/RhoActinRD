function ParametersInPrior = ParameterSetsInPrior(ParInds,FixedParSet)
    if (FixedParSet>0)
        nPerDim=200;
        nV=2;
    else
        nPerDim=60;
        nV=4;
    end
    Allps = zeros(nPerDim,nV);
    TotalNum = nPerDim^nV;
    if (FixedParSet==0)
        PBounds = [1 40; 0.1 10; 0 8; 0 100];
        for k=1:size(PBounds,1)
            Allps(:,k) = PBounds(k,1)+(1:nPerDim)/nPerDim*(PBounds(k,2)-PBounds(k,1));
        end
        [x1,x2,x3,x4]=ndgrid(Allps(:,1),Allps(:,2),Allps(:,3),Allps(:,4));
        EveryParameter = [0.7*ones(TotalNum,1) 0.4*ones(TotalNum,1) x1(:) ones(TotalNum,2) ...
            x2(:) x3(:) x4(:)];
    elseif (FixedParSet==1)
        PBounds = [1 40; 0.1 10];
        for k=1:size(PBounds,1)
            Allps(:,k) = PBounds(k,1)+(1:nPerDim)/nPerDim*(PBounds(k,2)-PBounds(k,1));
        end
        [x1,x2]=ndgrid(Allps(:,1),Allps(:,2));
        EveryParameter = [0.7*ones(TotalNum,1) 0.4*ones(TotalNum,1) x1(:) ones(TotalNum,2) ...
            x2(:) 1.1*ones(TotalNum,1) 50*ones(TotalNum,1)];
    elseif (FixedParSet==2)
        PBounds = [0 4; 0 100];
        for k=1:size(PBounds,1)
            Allps(:,k) = PBounds(k,1)+(1:nPerDim)/nPerDim*(PBounds(k,2)-PBounds(k,1));
        end
        [x3,x4]=ndgrid(Allps(:,1),Allps(:,2));
        EveryParameter = [[0.7 0.4 30 1 1 1].*ones(TotalNum,6) x3(:) x4(:)];
    end
    MonomerClock = EveryParameter(:,6)./EveryParameter(:,4)+...
        EveryParameter(:,6)./EveryParameter(:,5)+EveryParameter(:,3);
    NucRate_B = EveryParameter(:,7)./(MonomerClock.*EveryParameter(:,6));
    NucRate_Rho = EveryParameter(:,8)./(MonomerClock.*EveryParameter(:,6));
    EveryParameter(:,7)=NucRate_B;
    EveryParameter(:,8)=NucRate_Rho;
    EveryParameter=AddParams(EveryParameter);
    
    load('TC_uStimInducSustain.mat','TC_Sustainable')
    InPrior = TC_Sustainable.predictFcn(EveryParameter(:,ParInds));
    ParametersInPrior = EveryParameter(InPrior,:);
end