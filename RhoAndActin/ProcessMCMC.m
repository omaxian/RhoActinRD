Names=["Starfish" "nmy" "nmy-pfn" "nmy-cyk"];
Symbols = ["o" "s" ">" "d"];
for jName=1:4
EmType=Names(jName);
load(strcat(EmType,'MCMCRunBNEq.mat'))
nP=nParams;
% For MCMC
StartIndex=nSamp/2+1;
EndIndex=nSamp;
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
%AllParameters(9,:)=AllParameters(3,:)+AllParameters(6,:)./AllParameters(4,:)+...
%    AllParameters(6,:)./AllParameters(5,:);
%AllParameters(10,:)=AllMeanActins;

LogLikelihood=AllDiffNorms+AllExSizeErs;
MaxDiff = max(AllDiffNorms);
% re-order indices
[v,x]=sort(LogLikelihood,'descend');
inds=x;
%inds(~Accepted(inds))=[];
inds=inds(end-99:end);
PL1=AllParameters(7,inds)<0.03;
% AllParameters(3,:)=AllParameters(3,:)+...
%     2*AllParameters(6,:)./AllParameters(4,:);
TotalLifetime = AllParameters(6,:)./AllParameters(4,:)+...
    AllParameters(6,:)./AllParameters(5,:)+AllParameters(3,:);
NFilBulk = AllParameters(7,:).*TotalLifetime;
MeanBulk = NFilBulk.*AllParameters(6,:);
NFilRho = AllParameters(8,:).*TotalLifetime;
MeanRho = NFilRho.*AllParameters(6,:);
ChaseSpeed = AllParameters(8,:).*AllParameters(6,:);
% Estimate mesh size 
load('MeshSizeEstimateData.mat')
Ls=AllParameters(6,:);
MeshFactor=0*Ls;
for j=1:length(Ls)
    Lfil=Ls(j);
    if (Lfil < Lens(1))
        MeshFactor(j)=fLens(1);
    elseif (Lfil>Lens(end))
        MeshFactor(j)=fLens(end);
    else
    [~,ind1]=min(abs(Lens-Lfil));
    L1 = Lens(ind1);
    if (L1 < Lfil)
        ind2 =ind1+1;
    else
        ind2=ind1-1;
    end
    % Linear interpolation
    MeshFactor(j) = fLens(ind1)+(fLens(ind2)-fLens(ind1))./...
        (Lens(ind2)-Lens(ind1)).*(Lfil-Lens(ind1));
    end
end

L=20;
Hetero = MeshFactor.*L^2./(MeanBulk+MeanRho);
AllParameters(9,:)=MeanBulk;
AllParameters(10,:)=MeanRho;
AllParameters(11,:)=Hetero;
xIndex = [1 3 4 7 9];
yIndex = [2 6 6 8 10];
xLabels = ["$k_\textrm{off}^{(0)}$" "$T_\textrm{fil}$" "$\nu_\textrm{p}$" ...
    "$q_b$" "$T_p$"];
yLabels = ["$k_\textrm{inh}$" "$\ell_\textrm{max}$" "$\ell_\textrm{max}$" ...
    "$q_\rho$" "$T_\textrm{fil}$"];
xLimits = [0.25 1.22 0 30 0 5 0 0.5 0 2.5];
yLimits = [0.2 1.5 0 20 0 20 0 1.5 0 7.5];
%tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
if (jName <4)
scatter3(AllParameters(9,inds),AllParameters(10,inds),...
    AllParameters(11,inds),Symbols(jName),'filled')
hold on
else
scatter3(AllParameters(9,inds(~PL1)),AllParameters(10,inds(~PL1)),...
    AllParameters(11,inds(~PL1)),Symbols(jName),'filled')
set(gca,'ColorOrderIndex',jName)
scatter3(AllParameters(9,inds(PL1)),AllParameters(10,inds(PL1)),...
    AllParameters(11,inds(PL1)),Symbols(jName))
end
for iP=2:-1
%nexttile
subplot(1,3,iP-1)
%scatter(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds),...
%    20,LogLikelihood(inds),Symbols(jName),'filled')
if (jName <4 )
scatter(AllParameters(xIndex(iP),inds),AllParameters(yIndex(iP),inds),...
    Symbols(jName),'filled')
hold on
else
scatter(AllParameters(xIndex(iP),inds(PL1)),AllParameters(yIndex(iP),inds(PL1)),...
    Symbols(jName))
hold on
set(gca,'ColorOrderIndex',jName)
scatter(AllParameters(xIndex(iP),inds(~PL1)),AllParameters(yIndex(iP),inds(~PL1)),...
    Symbols(jName),'filled')
end

pbaspect([1 1 1])
clim([min(LogLikelihood) 1])
colormap(flipud(turbo))
% box on
xlabel(xLabels(iP))
ylabel(yLabels(iP))
xlim([xLimits(2*iP-1) xLimits(2*iP)])
ylim([yLimits(2*iP-1) yLimits(2*iP)])
box on
end

if (0)
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
if (jName==5)
    CandInds=CandInds(Nnzs==numNonZero & AllParameters(7,:)'<0.03);
else
    CandInds=CandInds(Nnzs==numNonZero);
end
[~,ind]=min(abs(LogLikelihood(CandInds)));
ind=CandInds(ind);


%close all;
ps=AllParameters(1:end,ind);
ExSizesAll=[];
nNz=0;
TotActin=0;
for seed=1:nSeed
    tic
    Stats=RhoAndActinBasalNuc(ps,seed,nNz==0);
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
imagesc(Stats.rSim,Stats.tSim,XCorAvg)
xlim([0 max(Uvals)])
ylim([-120 120])
clim([-1 1])
colorbar
colormap turbo
xlabel('$\Delta r$ ($\mu$m)')

if (jName>1)
nexttile
plot(dsHist/2:dsHist:400,xp)
hold on
plot(dsHist/2:dsHist:400,SizeHist)
xlabel('Ex size ($\mu$m$^2$)')
ylabel('pdf')
xlim([0 200])
end
end
end