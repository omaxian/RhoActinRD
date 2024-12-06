% Assumes that actin polymerizes at some rate, stays at a fixed length for
% some rate, and then shrinks at some rate
function [XfEnd,nPerFil,TimeAtMax,AddedPts,RemovedPts] = UpdateActin(Xf,nPerFil,...
    systime,TimeAtMax,FullLifetime,GrowAmt,ShrinkAmt,MaxLength,ds)
    nFil=length(nPerFil);
    StartIndex = [0;cumsum(nPerFil)];
    nTot=StartIndex(end);
    RemovePoints = zeros(nTot,1);
    RemoveFibs = zeros(nFil,1);
    nMaxMon=floor(MaxLength/ds+1e-5)+1;
    XfEnd=nan*zeros(nMaxMon*nFil,2);
    AddedPts = nan*zeros(nMaxMon*nFil,2);
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
            % If shrink amount is a non-integer, randomly round up or down
            thisShrink = floor(ShrinkAmt)+1.0*(rand < mod(ShrinkAmt,1));
            nToRemove = min(thisShrink,nPerFil(iF));
            RemovePoints(stMon+1:stMon+nToRemove)=1;
            nPerFil(iF)=nPerFil(iF)-nToRemove;
            XfEnd((iF-1)*nMaxMon+1:(iF-1)*nMaxMon+nPerFil(iF),:)...
                =XThis(nToRemove+1:end,:);
            if (nPerFil(iF)==0)
                RemoveFibs(iF)=1;
            end
        elseif (TimeAtMax(iF)==inf)
            % Polymerizing
            % If grow amount is a non-integer, randomly round up or down
            thisGrow = floor(GrowAmt)+1.0*(rand < mod(GrowAmt,1));
            nToAdd = min(nMaxMon-nPerFil(iF),thisGrow);
            if (nPerFil(iF)>1)
                tau=XThis(2,:)-XThis(1,:);
            else 
                tau=randn(1,2);
                tau=tau/norm(tau)*ds;
            end
            nPerFil(iF)=nPerFil(iF)+nToAdd;
            XfEnd((iF-1)*nMaxMon+1:(iF-1)*nMaxMon+nPerFil(iF),:)...
                =[XThis; (1:nToAdd)'*tau+XThis(end,:)];
            AddedPts((iF-1)*nMaxMon+1:(iF-1)*nMaxMon+nToAdd,:)...
                =(1:nToAdd)'*tau+XThis(end,:);
            if (nPerFil(iF)==nMaxMon)
                TimeAtMax(iF)=systime;
            end
        else
            XfEnd((iF-1)*nMaxMon+1:(iF-1)*nMaxMon+nPerFil(iF),:)=XThis;
        end
    end
    XfEnd=rmmissing(XfEnd);
    AddedPts=rmmissing(AddedPts);
    RemovedPts=Xf(RemovePoints==1,:);
    nPerFil(RemoveFibs==1)=[];
    TimeAtMax(RemoveFibs==1)=[];
end