% This file is used to sort through all the sample parameter sets
EveryParams=[];
EveryDiffNorm=[];
EveryExSizeEr=[];
EveryMeanActin=[];
EverynNz=[];
for seed=1:250
    try
    load(strcat('ScanFullICWtdk8_',num2str(seed),'.mat'))
    EveryParams=[EveryParams; AllParameters'];
    EveryDiffNorm=[EveryDiffNorm;AllDiffNorms];
    EveryExSizeEr=[EveryExSizeEr;AllExSizeErs];
    EveryMeanActin=[EveryMeanActin;AllMeanActins];
    EverynNz = [EverynNz;nNzs];
    catch
        seed
    end
end
inds1=find(EverynNz==2);
EveryParams=EveryParams(inds1,:);
EveryDiffNorm=EveryDiffNorm(inds1,:);
EveryExSizeEr=EveryExSizeEr(inds1,:);
EveryMeanActin=EveryMeanActin(inds1);
EverynNz=EverynNz(inds1);
[EveryParams,indrs]=uniquetol(EveryParams,'ByRows',true);
EveryDiffNorm=EveryDiffNorm(indrs,:);
EveryExSizeEr=EveryExSizeEr(indrs,:);
EveryMeanActin=EveryMeanActin(indrs);
EverynNz=EverynNz(indrs);
MonomerClock = EveryParams(:,6)./EveryParams(:,4)+EveryParams(:,3);
EveryParams(:,9)=MonomerClock;
TotalLifetime = EveryParams(:,6)./EveryParams(:,4)+...
    EveryParams(:,6)./EveryParams(:,5)+EveryParams(:,3);
NFilBulk = EveryParams(:,7).*TotalLifetime;
MeanBulk = NFilBulk.*EveryParams(:,6);
EveryParams(:,10)=MeanBulk;
NFilRho = EveryParams(:,8).*TotalLifetime;
MeanRho = NFilRho.*EveryParams(:,6);
EveryParams(:,11)=MeanRho;
TEr = EveryDiffNorm;
for iType=1:4
if (iType>1)
    TEr(:,iType)=EveryDiffNorm(:,iType)+EveryExSizeEr(:,iType);
    % % Penalize actin
    if (iType==2)
        ActinEr = (EveryMeanActin/0.5-1).^2.*(EveryMeanActin<0.5);
    elseif (iType==3)
        ActinEr = 2*(MonomerClock/12-1).^2.*(MonomerClock<12) + ...
            2*(EveryMeanActin/0.3-1).^2.*(EveryMeanActin> 0.3);
    elseif (iType==4)
        ActinEr = 4*(EveryMeanActin/0.3-1).^2.*(EveryMeanActin>0.3); 
    end
    TEr(:,iType) = TEr(:,iType) + ActinEr;
end
[vv,inds]=sort(TEr(:,iType),'descend');
inds=inds(end-50:end);
indreps(iType)=inds(end);
xIndex = [1 9 10 4 7];
yIndex = [2 6 11 6 8];
%xIndex = [1 3 4 7];
%yIndex = [2 6 6 8];
for iP=2:3
%nexttile
subplot(1,3,iP-1)
pltc{iP}(iType)=scatter(EveryParams(inds,xIndex(iP)),...
    EveryParams(inds,yIndex(iP)),20,'filled');
hold on
box on
uistack(pltc{iP}(iType),'bottom')
end
subplot(1,3,3)
counts=histcounts(EveryMeanActin(inds),0:0.1:0.9);
plot(0.05:0.1:0.9,counts,'-o')
hold on
end
% 
close all;
figure(2)
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact')
figure(3)
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact')
figure(4)
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact')
for iType=4
ps=EveryParams(indreps(iType),:);
ExSizesAll=[];
nNz=0;
TotActin=0;
for seed=1:nSeed
    tic
    Stats=RhoAndActinBasalNuc(ps,seed,nNz==0); % use nNz==1 for Fig 2 in paper
    toc
    % Compute the norm relative to the experiment and the
    % difference in the excitation size (for C. elegans only)
    if (Stats.XCor(1)==0)
    else
        ExSizesAll=[ExSizesAll;Stats.ExSizes];
    end
    % Cross correlation difference
    XCorEr = 1;
    if (Stats.XCor(1)~=0) 
        nNz=nNz+1;
        if (nNz==1)
            XCorAvg=Stats.XCor;
        else
            XCorAvg=XCorAvg+Stats.XCor;
        end
        TotActin=TotActin+Stats.MeanActin;
        tSimulated=Stats.tSim;
        rSimulated=Stats.rSim;
        if (nNz==numNonZero)
            break
        end
    end
end
% Compute errors 
XCorAvg=XCorAvg/nNz;
InterpolatedSim=ResampleXCor(XCorAvg,tSimulated,rSimulated,...
                    AllUvals{iType},Alldtvals{iType},...
                    max(AllUvals{iType})+1e-3,...
                    max(Alldtvals{iType})+1e-3);
XCorNorm=AllXCorsExp{iType}.^2;
ZeroEr = round(sum(XCorNorm(:)),1);
XCorEr = (InterpolatedSim-AllXCorsExp{iType}).^2;
XCorEr = sum(XCorEr(:))/ZeroEr
if (iType>1)
    xp=histcounts(ExSizesAll,0:dsHist:400);
    xp=xp/(sum(xp)*AlldsHist{iType});
    Wts = (dsHist/2:dsHist:400);
    ExSizeDiff = sum((xp-AllSizeHist{iType})...
        .*(xp-AllSizeHist{iType}).*Wts)...
        /sum(AllSizeHist{iType}.*AllSizeHist{iType}.*Wts) %L^2 norm
else
    ExSizeDiff=0;
end
MActin=TotActin/nNz

figure(5)
if (iType==2)
 tiledlayout(1,6,'Padding', 'none', 'TileSpacing', 'compact')
end
nexttile
imagesc(rSimulated,tSimulated,XCorAvg)
xlim([0 max(Uvals)])
ylim([-120 120])
clim([-1 1])
xlabel('$\Delta r$ ($\mu$m)')

nexttile
imagesc(AllUvals{iType},Alldtvals{iType},AllXCorsExp{iType})
xlim([0 max(Uvals)])
ylim([-120 120])
clim([-1 1])
colormap turbo
xlabel('$\Delta r$ ($\mu$m)')
yticklabels('')

%nexttile
figure(6)
if (iType>1)
plot(dsHist/2:dsHist:400,xp)
hold on
set(gca,'ColorOrderIndex',iType-1)
plot(dsHist/2:dsHist:400,AllSizeHist{iType},':')
xlabel('Ex size ($\mu$m$^2$)')
ylabel('pdf')
xlim([0 200])
end
end

