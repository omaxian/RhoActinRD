% Makes Fig. 4 - identifying parameter sets where the error
% in the cross correlation is the same but the log likelihood is 
% different
% Test set where we know the true result of fwd model
rng(0);
nP=5;
if (nP==2)
    ParInds = [3 6];
elseif (nP==5)
    ParInds = [2 3 6 10 11];
end
load(strcat('UniformResults/Scan',num2str(nP),'NewOneInPrior.mat'))
load('PrincipalComponentsXCor.mat')
load('AllRepeats5.mat')
IndexSim=4;
pTrue=AddParams(Allps(IndexSim,:));
TruePar = pTrue(ParInds);
TestXCorsAll = XCorsPars(:,:,IndexSim);
nLay=[20 40 20];
nW=1;
Noise = 0;
nModes=5;
RegZations=1e-4;

nTest=5000;
testinds = randperm(length(AllParameters),nTest);
Ers=zeros(1,nTest);
nSampsTrue = size(TestXCorsAll,2);
for p=1:nSampsTrue
    D = AllXCors(:,testinds)-TestXCorsAll(:,p);
    Ers = Ers+sum(D.*D)/nSampsTrue;
end
Ers = sqrt(Ers./sum(mean(TestXCorsAll,2).*mean(TestXCorsAll,2)));
SamplesXCors=U'*AllXCors(:,testinds);
SamplesXCors = SamplesXCors(1:nModes,:);
pChk = AllParameters(ParInds,testinds);

% Sort the errors and take some evenly spaced ones
InputVec = [SamplesXCors'.*ones(nTest,1) repmat(TruePar,nTest,1)];
LikelihoodToEvidence=zeros(nTest,1);
nSeed=2;
for seed=1:nSeed
load(strcat('NewTC',num2str(seed-1),'_XCorNoise',num2str(Noise),...
    'PCA',num2str(nModes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
[yfit,scores] = trainedClassifier.predictFcn(InputVec); 
LikelihoodToEvidence = LikelihoodToEvidence+1/nSeed*scores(:,2)./(1-scores(:,2));
end
%nexttile
scatter(Ers,log10(LikelihoodToEvidence),20,'filled');
%tc = log10(LikelihoodToEvidence) > -10;
%r = corrcoef(Ers(tc),log10(LikelihoodToEvidence(tc)));
%legend(strcat('$r=$',num2str(round(r(1,2)*100)/100)))
xlabel('$||x-\bar x(\theta)||$','interpreter','latex')
ylabel('$\log_{10} \langle \hat r(x,\theta) \rangle$','interpreter','latex')
%ylim([-10 5])

% Find the first xcor error where there are low and high probability
% observations
ds=0.1;
st = min(Ers)-ds;
while st < 1
    st=st+ds/5;
    ErRange = [st st+ds];
    Inds = find(Ers > ErRange(1)& Ers<ErRange(2));
    LikRange = max(log10(LikelihoodToEvidence(Inds)))...
        -min(log10(LikelihoodToEvidence(Inds)));
    if (~isempty(Inds)&&LikRange>3)
        break
    end
end
hold on
plot((ErRange(1))*[1 1],ylim,':k')
plot((ErRange(2))*[1 1],ylim,':k')

% Find the lik where there are low and high XCor errors
ds=1;
st = max(log10(LikelihoodToEvidence))-ds;
while st > 0
    st=st-ds/10;
    LikRange = [st st+ds];
    IndsPlus = find(log10(LikelihoodToEvidence) > LikRange(1)& ...
        log10(LikelihoodToEvidence)<LikRange(2));
    ErRange = max(Ers(IndsPlus))-min(Ers(IndsPlus));
    if (ErRange>0.5)
        break
    end
end

%ErRange = [1.065 1.085];
plot(xlim,(LikRange(1))*[1 1],':k')
plot(xlim,(LikRange(2))*[1 1],':k')
for p=1:2
figure(1+p)
if (p==1)
% Fixed L
[~,orderplot]=sort(Ers(IndsPlus),'ascend');
%orderplot=orderplot(round(linspace(1, length(orderplot), 5)));
orderplot = orderplot([1 12 26 34 35]);
IndsPlot=IndsPlus(orderplot);
tiledlayout(1,length(IndsPlot)+1,'Padding', 'none', ...
    'TileSpacing', 'compact')
else
[~,orderplot]=sort(LikelihoodToEvidence(Inds),'descend');
%orderplot=orderplot([1 3 5 7 8]);
%orderplot=orderplot(round(linspace(1, length(orderplot), 5)));
orderplot = orderplot([1 9 13 19 20]);

IndsPlot = Inds(orderplot);
% Fixed er
tiledlayout(length(IndsPlot),1,'Padding', 'none', ...
    'TileSpacing', 'compact') 
end
if (p==1)
nexttile
imagesc(0:0.5:10,-120:5:120,...
    reshape(mean(TestXCorsAll,2),121,21))
clim([-1 1])
colormap turbo
pbaspect([1 1 1])
xticklabels('')
yticklabels('')
title('Mean XCor')
end
for k=1:length(IndsPlot)
nexttile
imagesc(0:0.5:10,-120:5:120,...
    reshape(AllXCors(:,testinds(IndsPlot(k))),121,21))
clim([-1 1])
colormap turbo
if (p==1)
title(strcat('$\bar e=',num2str(round(Ers(IndsPlot(k)),1)),'$'),'interpreter','latex')
else
ylabel(strcat('$\bar \ell=',num2str(round(log10(LikelihoodToEvidence(IndsPlot(k))),1)),'$'),'interpreter','latex')
end    
pbaspect([1 1 1])
xticklabels('')
yticklabels('')
end
figure(4)
if (p==1)
tiledlayout(2,1,'Padding', 'none', 'TileSpacing', 'compact');
end
nexttile
CRng=U'*AllXCors(:,testinds(IndsPlot));
TestXCorsHat=U'*TestXCorsAll;
TestXCorsHat=TestXCorsHat(1:nModes,:);
for c=1:nModes
if (p==2)
scatter(fliplr(CRng(c,:)),c*ones(1,length(IndsPlot)),36,...
    flipud(log10(LikelihoodToEvidence(IndsPlot))),'filled')
else
scatter(fliplr(CRng(c,:)),c*ones(1,length(IndsPlot)),36,...
    fliplr(Ers(IndsPlot))','filled')
end
hold on
ctrBin=round(mean(TestXCorsHat(c,:)));
dx=0.5;
bedge=ctrBin-5:dx:ctrBin+5;
xpl = 1/2*(bedge(1:end-1)+bedge(2:end));
True(c,:) =histcounts(TestXCorsHat(c,:),bedge);
True(c,:) =True(c,:)/(dx*size(TestXCorsHat,2));
plot(xpl,True(c,:)+c,'-k')
end
turbo(100)
c=ans(30:85,:);
colormap(gca,c)
if (p==1)
    colormap(gca,flipud(c))
end
colorbar
box on
pbaspect([1 1 1])
xlabel('$\gamma_j$','interpreter','latex')
ylabel('$j$','interpreter','latex')
ylim([1 6])
yticks(1:2:5)
end
