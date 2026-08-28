% Experimental data plot in PCA and parameter space
for offset=0:5:15
nP=5;
nPv=5;
nLay=[20 40 20];
nW=1;
%offset=0;
nModesPlot=10;
nPerPlot = 5;
nTrial = 4;
Noise = 0;
if (offset>=10)
    nModes=10;
else
    nModes=5;
end
if (mod(offset,10)==0)
    RegZations=1e-4;
else
    RegZations=1e-3;
end
nToCheck=20000;
if (offset==0)
figure;
tiledlayout(4,5,'Padding', 'none', 'TileSpacing', 'compact');
end
rng(0);
load('XCor_Worm.mat')
load('PrincipalComponentsXCor.mat')
if (nP==2)
    ParInds = [3 6];
elseif (nP==5)
    ParInds = [2 3 6 10 11];
end
TestXCor=DataXCor(:);
TestXCorTr=U'*TestXCor;

% Load data and classifier
for iTri=1:nTrial
if (iTri<=nTrial/2)
load(strcat('NewTC0_XCorNoise',num2str(Noise),...
    'PCA',num2str(nModes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
else
load(strcat('NewTC1_XCorNoise',num2str(Noise),...
    'PCA',num2str(nModes),...
    'Lam',num2str(RegZations),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP',num2str(nP),'.mat'),...
    'trainedClassifier')
end
if (iTri==1)
load(strcat('UniformResults/Scan',num2str(nPv),'NewOneInPrior.mat'))
end
% Choose samples
inds=randperm(size(AllParameters,2),nToCheck);
ScanParams=AllParameters(:,inds);
ScanXCors=AllXCors(:,inds);
% ABC
er=DataXCor(:)-AllXCors;
sqer=sum(er.*er);
[~,minind]=min(sqer);
XCorMin = AllXCors(:,minind);
ParamsMin=AllParameters(:,minind);

TestCors=U'*ScanXCors;
XCorMinHat = U'*XCorMin;
% Parameters with data fixed
InputVec_Par = [repmat(TestXCorTr(1:nModes)',nToCheck,1) ScanParams(ParInds,:)'];
[~,scores] = trainedClassifier.predictFcn(InputVec_Par); 
LikelihoodToEvidence = scores(:,2)./(1-scores(:,2));
[~,plorder]=sort(LikelihoodToEvidence);
if (iTri==1)
if (nP==2)
    scatter(ScanParams(ParInds(1),:),ScanParams(ParInds(2),:),12,log10(LikelihoodToEvidence),'filled')
else
    % nexttile(offset+1)
    % scatter(AllParameters(ParInds(1),plorder),AllParameters(ParInds(5),plorder),...
    %     12,log10(LikelihoodToEvidence(plorder)),'filled')
    % g=clim;
    % clim([-3 ceil(max(g))])
    % t=turbo(100);
    % c=t(30:85,:);
    % colormap(gca,c)
    % title('Parameter space')
    % xlabel('$k_\textrm{GAP}$','interpreter','latex')
    % ylabel('$f_\rho$','interpreter','latex')
    % pbaspect([1 1 1])
    % nexttile(offset+2)
    % scatter(AllParameters(ParInds(4),plorder),AllParameters(ParInds(5),plorder),...
    %     12,log10(LikelihoodToEvidence(plorder)),'filled')
    % hold on
    % g=clim;
    % clim([-3 ceil(max(g))])
    % t=turbo(100);
    % c=t(30:85,:);
    % colormap(gca,c)
    % title('Parameter space')
    % xlabel('$f_b$','interpreter','latex')
    % ylabel('$f_\rho$','interpreter','latex')
    %pbaspect([1 1 1])
    nexttile(offset+1)
    scatter(ScanParams(ParInds(2),plorder),ScanParams(ParInds(3),plorder),...
        12,log10(LikelihoodToEvidence(plorder)),'filled')
    hold on
    scatter(ParamsMin(ParInds(2)),ParamsMin(ParInds(3)),'ko','filled')
end
g=clim;
clim([-3 ceil(max(g))])
t=turbo(100);
c=t(30:85,:);
colormap(gca,c)
title('Parameter space')
xlabel('$T_\textrm{fil}$','interpreter','latex')
ylabel('$\ell_\textrm{max}$','interpreter','latex')
pbaspect([1 1 1])
end
if (sum(isinf(LikelihoodToEvidence))>0)
    LikelihoodToEvidence(isinf(LikelihoodToEvidence))=1e3;
    keyboard
end
if (iTri==1)
    [~,plorder]=sort(LikelihoodToEvidence);
    IndexMax = inds(plorder(end));
    if (mod(IndexMax,2)==0)
        IndexMaxp1=IndexMax-1;
    else
        IndexMaxp1=IndexMax+1;
    end
    if (norm(AllParameters(:,IndexMax)-AllParameters(:,IndexMaxp1)) > 1e-12)
        keyboard
    end
    XCorMaxPost = 1/2*(AllXCors(:,IndexMax)+AllXCors(:,IndexMaxp1));
    XCorMaxPostHat = U'*XCorMaxPost;
    nexttile(offset+2)
    scatter(TestCors(1,plorder),TestCors(2,plorder),12,...
        log10(LikelihoodToEvidence(plorder)),'filled')
    hold on
    scatter(TestXCorTr(1,:),TestXCorTr(2,:),50,'ko','filled')
    g=clim;
    clim([-3 ceil(max(g))])
    t=turbo(100);
    c=t(30:85,:);
    colormap(gca,c)
    title('PCA space')
    xlabel('$\gamma_1$','interpreter','latex')
    ylabel('$\gamma_2$','interpreter','latex')
    pbaspect([1 1 1])
end

% Histograms
for j=1:nModesPlot
ctrBin=round(mean(TestXCorTr(j,:)));
dx=0.5;
bedge=ctrBin-5:dx:ctrBin+8;
xpl = 1/2*(bedge(1:end-1)+bedge(2:end));
AllCounts = zeros(length(xpl),1);
BinNum = 1+floor((TestCors(j,:)-bedge(1))/(bedge(2)-bedge(1)));
BinNum(BinNum<1)=1;
BinNum(BinNum>length(xpl))=length(xpl);
% Put the weights into the bins
TotalWeight = zeros(1,length(xpl));
for k=1:nToCheck    
    TotalWeight(BinNum(k))=TotalWeight(BinNum(k))+LikelihoodToEvidence(k);
end
TotalWeight=TotalWeight/(dx*sum(TotalWeight));

AllWts(j,:,iTri)=TotalWeight;
end
end

if (1)
% Summary plot
plcut=-1e-16;
for j=1:nModesPlot
ctrBin=round(TestXCorTr(j));
dx=0.5;
bedge=ctrBin-5:dx:ctrBin+8;
xpl = 1/2*(bedge(1:end-1)+bedge(2:end));
nexttile(offset+3+floor((j-1)/nPerPlot))
if (floor((j-1)/nPerPlot)==0)
    ylim([1 6])
    yticks(1:2:5)
else
    ylim([6 11])
    yticks(6:2:10)
end
plot(TestXCorTr(j)*[1 1],[j j+1],'-k')
hold on
plot(XCorMinHat(j)*[1 1],[j j+1],':k')
plot(XCorMaxPostHat(j)*[1 1],[j j+1],'-.r')
set(gca,'Colororderindex',2)
MeanWt = mean(AllWts(j,:,:),3);
ErWt = 2*std(AllWts(j,:,:),[],3)/sqrt(nTrial);
ErWt=ErWt(MeanWt>plcut);
xpl=xpl(MeanWt>plcut);
MeanWt=MeanWt(MeanWt>plcut);
plot(xpl,MeanWt+j)
cplot=get(gca,'ColorOrder');
hold on
fill([xpl, fliplr(xpl)], [j+MeanWt-ErWt fliplr(j+MeanWt+ErWt)],cplot(2,:),...
    'FaceAlpha', 0.2, 'linestyle', 'none');
pbaspect([1 1 1])
if (j==nModesPlot)
legend('Data','ABC min','Max post','$p(x|\theta)$','$\hat p(x|\theta)$','interpreter','latex','numColumns',2)
end
xlabel('$\gamma_j$','interpreter','latex')
ylabel('$j$','interpreter','latex')
end
end
nexttile(offset+3+floor(nModesPlot/nPerPlot))
imagesc(0:0.5:10,-120:5:120,reshape(XCorMaxPost,121,21))
colormap(gca,'turbo')
clim([-1 1])
pbaspect([1 1 1])
xlabel('$\Delta r$ ($\mu$m)','interpreter','latex')
ylabel('$\Delta t$ (s)','interpreter','latex')
title('Maximum posterior')
% nexttile(offset+4+floor(nModesPlot/nPerPlot))
% imagesc(0:0.5:10,-120:5:120,DataXCor)
% colormap(gca,'turbo')
% clim([-1 1])
% pbaspect([1 1 1])
% xlabel('$\Delta r$ ($\mu$m)','interpreter','latex')
% ylabel('$\Delta t$ (s)','interpreter','latex')
% title('Data')
end