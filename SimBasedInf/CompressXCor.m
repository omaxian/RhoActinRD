% Compress cross correlations using autoencoders and PCA
load('AllXCors.mat')
nTrain=10000;
nTest=1000;


%% Train the autoencode on the cross correlations
LatentSizes = [2 5 10 20 50];
if (0)
AllErs = zeros(1,length(LatentSizes));
for iP=1:length(LatentSizes)
maxep = 100*LatentSizes(iP);
autoenc = trainAutoencoder(AllXCors(:,1:nTrain),LatentSizes(iP),...
    'ShowProgressWindow',true,'MaxEpochs',maxep,...
        'SparsityRegularization',1,...
        'SparsityProportion',0.2);
NewSet = AllXCors(:,nTrain+1:nTrain+nTest);
EncodedNew = autoenc.encode(NewSet);
NewSetEnc = autoenc.decode(EncodedNew);
Er = NewSet-NewSetEnc;
AllErs(iP) = sum(sqrt(sum(Er.*Er))./sqrt(sum(NewSet.*NewSet)))/nTest;
AllAutoEncs{iP}=autoenc;
AllNewSets{iP}=NewSetEnc;
save('XCorAutoEncodersFull.mat','AllAutoEncs','LatentSizes');
end
figure;
nPl=10;
tiledlayout(2,3,'Padding', 'none', 'TileSpacing', 'compact')
indpick =255;%randperm(nTest,nPl);
pbaspect([1 1 1])
title('True','Fontweight','normal')
for k=1
    nexttile
    imagesc(0:0.5:10,-120:5:120,reshape(NewSet(:,indpick(k)),121,21))
    clim([-1 1])
    for jP=1:5
    nexttile
    imagesc(0:0.5:10,-120:5:120,reshape(AllNewSets{jP}(:,indpick(k)),121,21))
    clim([-1 1])
    colormap turbo
    xlabel('$\Delta r$ ($\mu$m)')
    ylabel('$\Delta t$ (s)')
    colorbar
    pbaspect([1 1 1])
    title(strcat(num2str(LatentSizes(jP)),' modes'),'interpreter','tex','Fontweight','normal')
    end
end
end

maxModes=1000;
X = AllXCors(:,1:nTrain);
[U,S,V]=svd(X,'econ');
Components=U;
Scores = S*V';
% for jj=1:5
% NewMat=NewMat+U(:,jj)*S(jj,jj)*V(:,jj)';
% end
NewSet = AllXCors(:,nTrain+1:nTrain+nTest);
NewScores = U'*NewSet;
pVals = [1:50 51:5:100 100:20:1000];
AllErs = zeros(size(pVals));
for iP=1:length(pVals)
    nM=pVals(iP);
    Rec = U(:,1:nM)*NewScores(1:nM,:);
    Er = Rec-NewSet;
    AllErs(iP) = sum(sqrt(sum(Er.*Er))./sqrt(sum(NewSet.*NewSet)))/nTest;
end
loglog(pVals,AllErs)