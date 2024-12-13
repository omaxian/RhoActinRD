EmType='nmy-cyk';
load(strcat(EmType,'MCMCRunL_NoWt.mat'))
nP=nParams;
% For MCMC
figure
nSoFar=iSamp;
Accepted=Accepted(1:nSoFar,:);
Accepted=Accepted(:)';
AllParametersES=AllParameters(:,1:nSoFar);
AllParameters=[];
for p=1:nWalker
    AllParameters = [AllParameters AllParametersES((p-1)*nP+1:p*nP,:)];
end
AllExSizeErs = reshape(AllExSizeErs(1:nSoFar,:),nSoFar*nWalker,1);
AllDiffNorms = reshape(AllDiffNorms(1:nSoFar,:),nSoFar*nWalker,1);
Nnzs = reshape(Nnzs(1:nSoFar,:),nSoFar*nWalker,1);
AllMeanActins = reshape(AllMeanActins(1:nSoFar,:),nSoFar*nWalker,1);
TotalLife = AllParameters(3,:) + 2*AllParameters(5,:)./AllParameters(4,:);
TotalSpace = AllParameters(5,:).*AllParameters(6,:);
AllParameters(7,:)=TotalLife;
AllParameters(8,:)=TotalSpace;

LogLikelihood=AllDiffNorms+AllExSizeErs;
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
%inds=x(end-49:end);
inds=x;
inds(~Accepted(inds))=[];
%AllParameters(6,:)=TwoMeanActins;
xIndex = [1 3 4 7];
yIndex = [2 5 6 8];
xLabels = ["$k_\textrm{off}^{(0)}$" "$T_\textrm{fil}$" "$\nu_p$" ...
    "$T_\textrm{fil}+2\ell_\textrm{max}/\nu_p$"];
yLabels = ["$k_\textrm{inh}$" "$\ell_\textrm{max}$" "$\bar q_0$" ...
    "$\bar{q}_0 \ell_\textrm{max}$"];
xLimits = [0.4 1.22 0 30 0 5 0 30];
yLimits = [0.2 1.1 0 20 0 1 0 8];
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
for iP=1:3
nexttile
%subplot(1,3,iP-1)
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
%scatter(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds),...
%    's','filled')
%hold on
% if (iP==4)
%     k=convhull(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds));
%     set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1)
%     plot(AllParameters(xIndex(iP),inds(k)),AllParameters(yIndex(iP),inds(k)));
% end
hold on
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
% figure;
% tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
% nexttile
% scatter(AllDifferencesModelExp(inds),TwoMeanActins(inds),'filled')
% xlabel('XCor Error')
% ylabel('Mean actin')
% % nexttile
% % scatter(AllDifferencesModelExp(inds),mean(AllNumEx(:,inds),'omitnan'),'filled')
% % box on
% % xlabel('XCor Error')
% % ylabel('Mean \# excitations')
% nexttile
% scatter(AllDifferencesModelExp(inds),mean(AllExSizeErs(:,inds),'omitnan'),'filled')
% xlabel('XCor Error')
% box on
% ylabel('Mean excitation size')
% 
% figure
% for iWalker=1:nWalker
% inds=(iWalker-1)*nSoFar+1:iWalker*nSoFar;
% plInds=inds(Accepted(inds)==1);
% Colors=1:length(plInds);
% scatter(AllDifferencesModelExp(plInds),TwoMeanSize(plInds),36,Colors,'filled')
% pause(1)
% end


figure
for iWalker=5:5:nWalker
inds=(iWalker-1)*nSoFar+1:iWalker*nSoFar;
plInds=inds(Accepted(inds)==1);
plot(mod(plInds-1,nSoFar),LogLikelihood(plInds),'-o')
hold on
end

figure
inds=1:nSoFar*nWalker;
plInds=inds(Accepted(inds)==1);
[~,order]=sort(LogLikelihood(plInds),'descend');
plInds=plInds(order);
Colors=mod(plInds,nSoFar);
scatter(AllDiffNorms(plInds),...
    AllExSizeErs(plInds),20,Colors,'filled')


[~,ind]=min(LogLikelihood-1.4);

ExSizesAll=[];
nNz=0;
TotActin=0;
for seed=1:nSeed
    tic
    Stats=RhoAndActin(AllParameters(:,ind),seed);
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
    end
end
% Compute errors 
XCorAvg=XCorAvg/nNz;
InterpolatedSim=ResampleXCor(XCorAvg,tSimulated,rSimulated,...
            Uvals,dtvals,max(Uvals)+1e-3,max(dtvals)+1e-3);
XCorEr = TotWts.*(InterpolatedSim-XCorsExp).^2;
XCorEr = sum(XCorEr(:))/ZeroEr
xp=histcounts(ExSizesAll,0:dsHist:400);
WtsEx=ones(1,length(xp));
xp=xp/(sum(xp)*dsHist);
ExSizeDiff = sum((xp-SizeHist).*(xp-SizeHist).*WtsEx)...
        /sum(SizeHist.*SizeHist.*WtsEx) %L^2 norm
MActin=TotActin/nNz

figure
imagesc(Uvals,dtvals,InterpolatedSim)
xlim([0 5])
ylim([-120 120])
clim([-1 1])
colorbar
colormap turbo
figure
imagesc(Uvals,dtvals,XCorsExp)
clim([-1 1])
colormap turbo
ylim([-120 120])
figure
plot(dsHist/2:dsHist:400,xp)
hold on
plot(dsHist/2:dsHist:400,SizeHist)