% Assumes that actin polymerizes at some rate, stays at a fixed length for
% some rate, and then shrinks at some rate
function [XfEnd,nPerFil,TimeAtMax,AddedPts,RemovedPts] = UpdateActin(Xf,nPerFil,...
    systime,TimeAtMax,FullLifetime,GrowAmt,ShrinkAmt,MaxLength,ds)
    nFil=length(nPerFil);
    StartIndex = [0;cumsum(nPerFil)];
    XfEnd=[];
    RemovedPts=[];
    AddedPts=[];
    RemoveInds=[];
    nMaxMon=MaxLength/ds+1;
    if (mod(nMaxMon,1)~=0)
        keyboard
    end
    for iF=1:nFil
        stMon = StartIndex(iF);
        % Identify the monomers
        XThis=Xf(stMon+1:stMon+nPerFil(iF),:);
        % Figure out what to do with the filament
        if (systime-TimeAtMax(iF)>FullLifetime)
            % Depolymerizing
            nToRemove = min(ShrinkAmt,nPerFil(iF));
            RemovedPts=[RemovedPts;XThis(1:nToRemove,:)];
            XThis(1:nToRemove,:)=[];
            nPerFil(iF)=nPerFil(iF)-nToRemove;
            if (nPerFil(iF)==0)
                RemoveInds=[RemoveInds;iF];
            end
        elseif (TimeAtMax(iF)==inf)
            % Polymerizing
            nToAdd = min(nMaxMon-nPerFil(iF),GrowAmt);
            try
                tau=XThis(2,:)-XThis(1,:);
            catch 
                tau=randn(1,2);
                tau=tau/norm(tau)*ds;
            end
            AddedPts=[AddedPts;(1:nToAdd)'*tau+XThis(end,:)];
            XThis=[XThis; (1:nToAdd)'*tau+XThis(end,:)];
            nPerFil(iF)=nPerFil(iF)+nToAdd;
            if (nPerFil(iF)==nMaxMon)
                TimeAtMax(iF)=systime;
            end
        end
        XfEnd=[XfEnd; XThis];
    end
    nPerFil(RemoveInds)=[];
    TimeAtMax(RemoveInds)=[];
end