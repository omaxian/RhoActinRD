EmType='nmy';
load(strcat(EmType,'MCMCRunPDE_Det.mat'))
nP=nParams;
%tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
% For MCMC
StartIndex=1;
EndIndex=iSamp;
nSoFar=EndIndex-StartIndex+1;
Accepted=Accepted(StartIndex:EndIndex,:);
Accepted=Accepted(:);
AllParametersES=AllParameters(:,StartIndex:EndIndex);
%SimSeeds=SimSeeds(:,StartIndex:EndIndex);
AllParameters=[];
AllSeeds=[];
for p=1:nWalker
    AllParameters = [AllParameters AllParametersES((p-1)*nP+1:p*nP,:)];
    %AllSeeds = [AllSeeds SimSeeds((p-1)*nSeed+1:p*nSeed,:)];
end
AllExSizeErs = reshape(AllExSizeErs(StartIndex:EndIndex,:),nSoFar*nWalker,1);
AllDiffNorms = reshape(AllDiffNorms(StartIndex:EndIndex,:),nSoFar*nWalker,1);
AllMeanActins = reshape(AllMeanActins(1:nSoFar,:),nSoFar*nWalker,1);
LogLikelihood=AllDiffNorms+AllExSizeErs;
 
inds=1:length(AllDiffNorms);
% Ap2Inds = find(AllParameters(6,:)>0.1);
% [~,inds2]=sort(AllDiffNorms(Ap2Inds));
% Ap2Inds=Ap2Inds(inds2);
% inds=setdiff(inds,Ap2Inds);
% re-order indices
[v,x]=sort(LogLikelihood(inds),'descend');
%inds=x(end-50:end);
inds=inds(x);
xIndex = [1 3 5 4];
yIndex = [2 4 6 7];
xLabels = ["$k_\textrm{off}^{(0)}$" "$q_b$" ...
    "$k_\textrm{off}^\textrm{(act)}$" "$q_\rho$"];
yLabels = ["$k_\textrm{inh}$" "$q_\rho$" ...
    "$D_a$" "$r$"];
xLimits = [PBounds(xIndex(1),:) PBounds(xIndex(2),:) ...
    PBounds(xIndex(3),:) PBounds(xIndex(4),:)];
yLimits = [PBounds(yIndex(1),:) PBounds(yIndex(2),:) ...
    PBounds(yIndex(3),:) PBounds(yIndex(4),:)];
tiledlayout(1,4,'Padding', 'none', 'TileSpacing', 'compact');
for iP=1:4
nexttile
% scatter(AllParameters(xIndex(iP),Ap2Inds),AllParameters(yIndex(iP),Ap2Inds),...
%      10,LogLikelihood(Ap2Inds),'>','filled')
% hold on
scatter(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds),...
    20,LogLikelihood(inds),'s','filled')
box on
hold on
clim([min(LogLikelihood) 2])
colormap(flipud(turbo))
xlabel(xLabels(iP))
ylabel(yLabels(iP))
xlim([xLimits(2*iP-1) xLimits(2*iP)])
ylim([yLimits(2*iP-1) yLimits(2*iP)])
end
%return
% 
% 
figure
for iWalker=5:5:nWalker
inds=(iWalker-1)*nSoFar+1:iWalker*nSoFar;
plInds=inds(Accepted(inds)==1);
plot(mod(plInds-1,nSoFar),LogLikelihood(plInds),'-o')
hold on
end
% 
figure
inds=1:nSoFar*nWalker;
plInds=inds(Accepted(inds)==1);
[~,order]=sort(LogLikelihood(plInds),'descend');
plInds=plInds(order);
Colors=mod(plInds,nSoFar);
scatter(AllDiffNorms(plInds),...
    AllExSizeErs(plInds),20,Colors,'filled')

CandInds=1:nSoFar*nWalker;
[vals,srtinds]=sort(LogLikelihood(CandInds));
CandInds=CandInds(srtinds);
%RecomputedLikelihood=zeros(nToCheck,1);
for jInd=1
ps=AllParameters(:,CandInds(jInd));
ExSizesAll=[];
TotActin=0;
for seed=1:nSeed
    if (Randomness)
        [Stats,st]=RhoAndActinPDEs_RandomNuc...
            (ps,dt,seed,1,Nuc0s,NucEns);
    else
        [Stats,st]=RhoAndActinPDEs(ps,dt,1);
    end
    % Compute the norm relative to the experiment and the
    % difference in the excitation size (for C. elegans only)
    if (Stats.XCor(1)==0)
    else
        ExSizesAll=[ExSizesAll;Stats.ExSizes];
    end
    % Cross correlation difference
    XCorEr = 1;
    if (Stats.XCor(1)~=0) 
        if (seed==1)
            XCorAvg=Stats.XCor;
        else
            XCorAvg=XCorAvg+Stats.XCor;
        end
        TotActin=TotActin+Stats.MeanActin;
        tSimulated=Stats.tSim;
        rSimulated=Stats.rSim;
    end
end
% Compute errors 
XCorAvg=XCorAvg/nSeed;
InterpolatedSim=ResampleXCor(XCorAvg,tSimulated,rSimulated,...
            Uvals,dtvals,max(Uvals)+1e-3,max(dtvals)+1e-3);
XCorEr = TotWts.*(InterpolatedSim-XCorsExp).^2;
XCorEr = sum(XCorEr(:))/ZeroEr
if (~(EmType=='Starfish'))
    xp=histcounts(ExSizesAll,0:dsHist:400);
    WtsEx=ones(1,length(xp));
    xp=xp/(sum(xp)*dsHist);
    ExSizeDiff = sum((xp-SizeHist).*(xp-SizeHist).*WtsEx)...
            /sum(SizeHist.*SizeHist.*WtsEx) %L^2 norm
else
    ExSizeDiff=0;
end
MActin=TotActin/nSeed;
end

figure
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
nexttile
imagesc(Uvals,dtvals,XCorsExp)
clim([-1 1])
colormap turbo
ylim([-120 120])

nexttile
imagesc(Uvals,dtvals,InterpolatedSim)
xlim([0 5])
ylim([-120 120])
clim([-1 1])
colorbar
colormap turbo

% nexttile
% plot(dsHist/2:dsHist:400,xp)
% hold on
% plot(dsHist/2:dsHist:400,SizeHist)