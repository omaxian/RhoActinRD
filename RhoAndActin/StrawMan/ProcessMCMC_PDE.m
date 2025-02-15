EmType='Starfish';
load(strcat(EmType,'MCMCRunPDE_All.mat'))
nP=nParams;
%tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
% For MCMC
StartIndex=1;
EndIndex=iSamp;
nSoFar=EndIndex-StartIndex+1;
Accepted=Accepted(StartIndex:EndIndex,:);
Accepted=Accepted(:)';
AllParametersES=AllParameters(:,StartIndex:EndIndex);
AllParameters=[];
for p=1:nWalker
    AllParameters = [AllParameters AllParametersES((p-1)*nP+1:p*nP,:)];
end
AllExSizeErs = reshape(AllExSizeErs(StartIndex:EndIndex,:),nSoFar*nWalker,1);
AllDiffNorms = reshape(AllDiffNorms(StartIndex:EndIndex,:),nSoFar*nWalker,1);
Nnzs = reshape(Nnzs(StartIndex:EndIndex,:),nSoFar*nWalker,1);
AllMeanActins = reshape(AllMeanActins(1:nSoFar,:),nSoFar*nWalker,1);

LogLikelihood=AllDiffNorms+AllExSizeErs;
%% 
MaxDiff = max(AllDiffNorms);
inds=1:length(AllDiffNorms);
Ap2Inds=find(abs(AllDiffNorms-1)<1e-10);
inds=setdiff(inds,Ap2Inds);
Ap3Inds=inds(AllDiffNorms(inds)>1+1e-5);
inds=setdiff(inds,Ap3Inds);
Ap4Inds=inds(AllParameters(1,inds)<0.55);
inds=setdiff(inds,Ap4Inds);
% re-order indices
[v,x]=sort(LogLikelihood,'descend');
inds=x;
inds(~Accepted(inds))=[];
%inds=inds(end-100:end);
%AllParameters(6,:)=TwoMeanActins;
xIndex = [1 3 5];
yIndex = [2 4 6];
xLabels = ["$k_\textrm{off}^{(0)}$" "$T_\textrm{fil}$" "$\bar{q}_{\rho}$"];
yLabels = ["$k_\textrm{inh}$" "$\ell_\textrm{max}$" "$\bar q_0$"];
xLimits = [0.2 1.25 0 0.5 0 1];
yLimits = [0.2 1.5 0 1.5 0 10];
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
for iP=1:3
nexttile
if (iP==1)
plot([0.55 0.55],[0 3],':k')
hold on
set(gca,'ColorOrderIndex',1)
end
% scatter(AllParameters(xIndex(iP),Ap2Inds),AllParameters(yIndex(iP),Ap2Inds),...
%     10,'d','filled','MarkerFaceColor',[0.83 0.82 0.78]); % gray
% hold on
% scatter(AllParameters(xIndex(iP),Ap3Inds),AllParameters(yIndex(iP),Ap3Inds),...
%     10,'^','filled', 'MarkerFaceColor',[1 0.6 0.78]) % pink
% scatter(AllParameters(xIndex(iP),Ap4Inds),AllParameters(yIndex(iP),Ap4Inds),...
%     10,LogLikelihood(Ap4Inds),'o','filled')
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
% figure
% for iWalker=5:5:nWalker
% inds=(iWalker-1)*nSoFar+1:iWalker*nSoFar;
% plInds=inds(Accepted(inds)==1);
% plot(mod(plInds-1,nSoFar),LogLikelihood(plInds),'-o')
% hold on
% end
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
CandInds=CandInds(Nnzs==numNonZero);
[~,ind]=min(abs(LogLikelihood(CandInds)));
ind=CandInds(ind);

ps=AllParameters(1:6,ind);
ExSizesAll=[];
nNz=0;
TotActin=0;
for seed=1:nSeed
    tic
    Stats=RhoAndActinPDEs(ps,seed,nNz==0);
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

nexttile
plot(dsHist/2:dsHist:400,xp)
hold on
plot(dsHist/2:dsHist:400,SizeHist)