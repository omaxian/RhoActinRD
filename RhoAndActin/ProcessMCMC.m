%load('GW25000.mat')
%load('FirstAttempt.mat')
load('MCMCActinCERun.mat')
nP=nParams;
% For MCMC
nSoFar=iSamp;
Accepted=Accepted(1:nSoFar,:);
Accepted=Accepted(:)';
AllParametersES=AllParameters(:,1:nSoFar);
AllParameters=[];
for p=1:nWalker
    AllParameters = [AllParameters AllParametersES((p-1)*nP+1:p*nP,:)];
end
AllAverageExSize = reshape(AllAverageExSize(:,1:nSoFar,:),nSeed,nSoFar*nWalker);
AllDiffNorms = reshape(AllDiffNorms(:,1:nSoFar,:),nSeed,nSoFar*nWalker);
AllMeanActins = reshape(AllMeanActins(:,1:nSoFar,:),nSeed,nSoFar*nWalker);
AllNumEx = reshape(AllNumEx(:,1:nSoFar,:),nSeed,nSoFar*nWalker);
TotalLife = AllParameters(3,:) + 2*AllParameters(5,:)./AllParameters(4,:);
TotalSpace = AllParameters(5,:).*AllParameters(6,:);
AllParameters(7,:)=TotalLife;
AllParameters(8,:)=TotalSpace;

AllDifferencesModelExp = mean(AllDiffNorms);
TwoMeanActins = mean(AllMeanActins);
TwoMeanSize = 2e-3*(mean(AllAverageExSize)-20).^2;
LogLikelihood=AllDifferencesModelExp+TwoMeanSize;
MaxDiff = max(AllDiffNorms);
inds=1:length(AllDifferencesModelExp);
Ap2Inds=find(abs(AllDifferencesModelExp-ZeroEr)<1e-10);
inds=setdiff(inds,Ap2Inds);
Ap3Inds=inds(AllDiffNorms(nSeed,inds)>ZeroEr+1e-5);
inds=setdiff(inds,Ap3Inds);
Ap4Inds=inds(AllParameters(1,inds)<0.55);
inds=setdiff(inds,Ap4Inds);
% re-order indices
[v,x]=sort(mean(AllDiffNorms(:,inds)),'descend');
%AllParameters(6,:)=TwoMeanActins;
xIndex = [1 3 4 7];
yIndex = [2 5 6 8];
xLabels = ["$k_\textrm{off}^{(0)}$" "$T_\textrm{fil}$" "$\nu_p$" ...
    "$T_\textrm{fil}+2\ell_\textrm{max}/\nu_p$"];
yLabels = ["$k_\textrm{inh}$" "$\ell_\textrm{max}$" "$\bar q_0$" ...
    "$\bar{q}_0 \ell_\textrm{max}$"];
xLimits = [0.4 1.22 0 30 0 5 0 30];
yLimits = [0.2 1.1 0 10 0 1 0 8];
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
for iP=2:4
nexttile
if (iP==1)
plot([0.55 0.55],[0 3],':k')
hold on
set(gca,'ColorOrderIndex',1)
end
scatter(AllParameters(xIndex(iP),Ap2Inds),AllParameters(yIndex(iP),Ap2Inds),...
    10,'d','filled','MarkerFaceColor',[0.83 0.82 0.78]); % gray
hold on
scatter(AllParameters(xIndex(iP),Ap3Inds),AllParameters(yIndex(iP),Ap3Inds),...
    10,'^','filled', 'MarkerFaceColor',[1 0.6 0.78]) % pink
scatter(AllParameters(xIndex(iP),Ap4Inds),AllParameters(yIndex(iP),Ap4Inds),...
    10,LogLikelihood(Ap4Inds),'o','filled')
scatter(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds),...
    20,LogLikelihood(inds),'s','filled')
colormap(flipud(turbo))
box on
xlabel(xLabels(iP))
ylabel(yLabels(iP))
xlim([xLimits(2*iP-1) xLimits(2*iP)])
ylim([yLimits(2*iP-1) yLimits(2*iP)])
if (iP==4)
    colorbar
end
end
%return;
figure;
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
scatter(AllDifferencesModelExp(inds),TwoMeanActins(inds),'filled')
xlabel('XCor Error')
ylabel('Mean actin')
nexttile
scatter(AllDifferencesModelExp(inds),mean(AllNumEx(:,inds),'omitnan'),'filled')
box on
xlabel('XCor Error')
ylabel('Mean \# excitations')
nexttile
scatter(AllDifferencesModelExp(inds),mean(AllAverageExSize(:,inds),'omitnan'),'filled')
xlabel('XCor Error')
box on
ylabel('Mean excitation size')