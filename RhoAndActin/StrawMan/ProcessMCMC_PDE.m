EmType='Starfish';
load(strcat(EmType,'MCMCRunPDE_SharpBox.mat'))
nP=nParams;
%tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
% For MCMC
StartIndex=1;
EndIndex=iSamp;
nSoFar=EndIndex-StartIndex+1;
Accepted=Accepted(StartIndex:EndIndex,:);
Accepted=Accepted(:);
AllParametersES=AllParameters(:,StartIndex:EndIndex);
SimSeeds=SimSeeds(:,StartIndex:EndIndex);
AllParameters=[];
AllSeeds=[];
for p=1:nWalker
    AllParameters = [AllParameters AllParametersES((p-1)*nP+1:p*nP,:)];
    AllSeeds = [AllSeeds SimSeeds((p-1)*nSeed+1:p*nSeed,:)];
end
AllExSizeErs = reshape(AllExSizeErs(StartIndex:EndIndex,:),nSoFar*nWalker,1);
AllDiffNorms = reshape(AllDiffNorms(StartIndex:EndIndex,:),nSoFar*nWalker,1);
Nnzs = reshape(Nnzs(StartIndex:EndIndex,:),nSoFar*nWalker,1);
AllMeanActins = reshape(AllMeanActins(1:nSoFar,:),nSoFar*nWalker,1);

AllnRts=zeros(nSoFar*nWalker,1);
AllnUnst=zeros(nSoFar*nWalker,1);
for j=1:nSoFar*nWalker
    [rts,stability] = PDERoots(AllParameters(:,j),0.1,20,100);
    AllnRts(j) = length(rts(:,1));
    AllnUnst(j) = sum(stability==-1);
end

LogLikelihood=AllDiffNorms+AllExSizeErs;
%% 
MaxDiff = max(AllDiffNorms);
inds=1:length(AllDiffNorms);
% Ap2Inds=find(AllnRts==3 & AllnUnst==2);
% Ap2Inds = find(AllParameters(6,:)>0.1);
% [~,inds2]=sort(AllDiffNorms(Ap2Inds));
% Ap2Inds=Ap2Inds(inds2);
% inds=setdiff(inds,Ap2Inds);
% Ap3Inds=find(AllnRts==3 & AllnUnst==1 & AllParameters(6,:)'<0.1);
% [~,inds3]=sort(AllDiffNorms(Ap3Inds));
% Ap3Inds=Ap3Inds(inds3);
% inds=setdiff(inds,Ap3Inds);
% Ap4Inds=find(AllnRts==1 & AllnUnst==1 & AllParameters(6,:)'<0.1);
% inds=setdiff(inds,Ap4Inds);
% [~,inds4]=sort(AllDiffNorms(Ap4Inds));
% Ap4Inds=Ap4Inds(inds4);
% re-order indices
[v,x]=sort(LogLikelihood(inds),'descend');
inds=inds(x);
% pChk=AllParameters(6,:)';
% inds(~Accepted(inds))=[];
% pChk(~Accepted(inds))=[];
% inds(pChk(inds)>0.1)=[];
% inds=inds(end-1000:end);
%AllParameters(6,:)=TwoMeanActins;
xIndex = [1 3 5];
yIndex = [2 4 6];
xLabels = ["$k_\textrm{off}^{(0)}$" "$\bar q_0$" ...
    "$k_\textrm{off}^{(act)}$"];
yLabels = ["$k_\textrm{inh}$" "$\bar{q}_{\rho}$" ...
    "$D^\textrm{(act)}$"];
xLimits = [PBounds(xIndex(1),:) PBounds(xIndex(2),:) ...
    PBounds(xIndex(3),:) ];
yLimits = [PBounds(yIndex(1),:) PBounds(yIndex(2),:) ...
    PBounds(yIndex(3),:) ];
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
for iP=1:3
nexttile
% scatter(AllParameters(xIndex(iP),Ap2Inds),AllParameters(yIndex(iP),Ap2Inds),...
%      10,LogLikelihood(Ap2Inds),'>','filled')
% hold on
% scatter(AllParameters(xIndex(iP),Ap3Inds),AllParameters(yIndex(iP),Ap3Inds),...
%      20,LogLikelihood(Ap3Inds),'^','filled')
% hold on
% scatter(AllParameters(xIndex(iP),Ap4Inds),AllParameters(yIndex(iP),Ap4Inds),...
%      10,LogLikelihood(Ap4Inds),'o','filled')
scatter(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds),...
    20,LogLikelihood(inds),'s','filled')
box on
%scatter(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds),...
%    'd','filled')
%box on
hold on
% if (iP==4)
%     k=convhull(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds));
%     set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1)
%     plot(AllParameters(xIndex(iP),inds(k)),AllParameters(yIndex(iP),inds(k)));
% end
% hold on
clim([min(LogLikelihood) 2])
colormap(flipud(turbo))
% box on
xlabel(xLabels(iP))
ylabel(yLabels(iP))
xlim([xLimits(2*iP-1) xLimits(2*iP)])
ylim([yLimits(2*iP-1) yLimits(2*iP)])
% if (iP==4)
%     colorbar
% end
end
% 
% figure
% for iWalker=1:nWalker
% inds=(iWalker-1)*nSoFar+1:iWalker*nSoFar;
% plInds=inds(Accepted(inds)==1);
% Colors=1:length(plInds);
% scatter(AllDifferencesModelExp(plInds),TwoMeanSize(plInds),36,Colors,'filled')
% pause(1)
% end
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
CandInds=CandInds(Nnzs==numNonZero  & AllParameters(6,:)'<0.1);
[~,ind]=min(abs(LogLikelihood(CandInds)));
ind=CandInds(ind);

ps=AllParameters(:,ind);
ExSizesAll=[];
nNz=0;
TotActin=0;
for seed=1:nSeed
    tic
    Stats=RhoAndActinPDEs(ps,dt,[],1);
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
MActin=TotActin/nNz

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