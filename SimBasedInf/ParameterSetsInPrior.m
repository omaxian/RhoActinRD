% Sample parameter sets from the prior
function [EveryParameter,InPrior] = ParameterSetsInPrior(ParInds)
    if (length(ParInds)==4)
        nPerDim=50;
        nV=4;
    elseif (length(ParInds)==5)
        nPerDim=25;
        nV = 5;
        if (nargin==2)
            nV = 3;
            nPerDim=80;
        end
    end
    Allps = zeros(nPerDim,nV);
    TotalNum = nPerDim^nV;
    if (length(ParInds)==4)
        PBounds = [1 40; 0.1 10; 0 8; 0 100];
        for k=1:size(PBounds,1)
            Allps(:,k) = PBounds(k,1)+(0:nPerDim-1)/(nPerDim-1)*(PBounds(k,2)-PBounds(k,1));
        end
        [x1,x2,x3,x4]=ndgrid(Allps(:,1),Allps(:,2),Allps(:,3),Allps(:,4));
        EveryParameter = [0.7*ones(TotalNum,1) 0.4*ones(TotalNum,1) x1(:) ones(TotalNum,2) ...
            x2(:) x3(:) x4(:)];
    elseif (length(ParInds)>4)
        PBounds = [0.2 0.6; 1 60; 0.1 10; 0 8; 0 150];
        for k=1:size(PBounds,1)
            Allps(:,k) = PBounds(k,1)+(0:nPerDim-1)/(nPerDim-1)*(PBounds(k,2)-PBounds(k,1));
        end
        [x1,x2,x3,x4,x5]=ndgrid(Allps(:,1),Allps(:,2),Allps(:,3),Allps(:,4),Allps(:,5));
        EveryParameter = [0.7*ones(length(x1(:)),1) x1(:) x2(:) ones(length(x1(:)),2) ...
            x3(:) x4(:) x5(:)];
        EveryParameter(EveryParameter(:,end-1)>1.6./EveryParameter(:,2),:) = [];
        EveryParameter(EveryParameter(:,end)> 150 - 140/0.4*(EveryParameter(:,2)-0.2),:) = [];
    end
    MonomerClock = EveryParameter(:,6)./EveryParameter(:,4)+...
        EveryParameter(:,6)./EveryParameter(:,5)+EveryParameter(:,3);
    NucRate_B = EveryParameter(:,7)./(MonomerClock.*EveryParameter(:,6));
    NucRate_Rho = EveryParameter(:,8)./(MonomerClock.*EveryParameter(:,6));
    EveryParameter(:,7)=NucRate_B;
    EveryParameter(:,8)=NucRate_Rho;
    EveryParameter=AddParams(EveryParameter);
    
    if (ParInds(1)~=2)
        load('TC_uStimInducSustain.mat','TC_Sustainable')
    else
        load('TC_rgInducSustain.mat','TC_Sustainable')
    end
    InPrior = TC_Sustainable.predictFcn(EveryParameter(:,ParInds));
    %ParametersInPrior = EveryParameter(InPrior,:);
end