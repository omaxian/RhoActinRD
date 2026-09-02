% This is the main function to train the classifier (it calls
% trainClassifier.m)
cd(fileparts(which(mfilename)));
addpath(genpath('..'))
for seed=1
rng(seed);
nP=5;
load(strcat('UniformResults/Scan',num2str(nP),'NewOneInPrior.mat'))
nSampAll=size(AllParameters,2);
inds=randperm(nSampAll,nSampAll);
AllParameters=AllParameters(:,inds);
AllXCors=AllXCors(:,inds);
load('PrincipalComponentsXCor.mat')
nTrial = size(AllParameters,2);
if (nP==2)
    ParInds=[3 6]; % parameters that vary
elseif (nP==4)
    ParInds=[3 6 10 11]; % parameters that vary
elseif (nP==5)
    ParInds=[2 3 6 10 11]; % parameters that vary
end
nNodes=[20 40 20];
nlev = [0]; % noise level
Lams = [0]; % L2 regularization strength
for iL=1:length(Lams)
for iW=[1]
SizesToTry = [2]; % number of PCA modes
for jEncSize=1:length(SizesToTry)
% Load data
XCorScores = U'*AllXCors; % These are the scores
MeanSqScore=mean(abs(XCorScores.^2),2);
numModes=SizesToTry(jEncSize);
XCorScores = XCorScores(1:numModes,:);
MeanSqScore=MeanSqScore(1:numModes)';
Noise = nlev(iL)*SingVals(1:numModes).*randn(numModes,nTrial);
XCorScores = XCorScores+Noise;

nPVary=length(ParInds);
nWrongPer=iW;
ClassifData_XC=zeros((1+nWrongPer)*nTrial,nPVary+1+size(XCorScores,1));
for j=1:nTrial
    theta0 = AllParameters(ParInds,j)';
    x0 = XCorScores(:,j)';
    % Check distance in PCA space
    PCADist = x0-XCorScores';
    PCADist = sum((PCADist.^2)./MeanSqScore,2);
    % How many nearby points are there
    nNearby(j)=sum(PCADist<0.01);
    ClassifData_XC((nWrongPer+1)*j-nWrongPer,:) = [x0 theta0 1];
    % Some other parameter not theta0 and not in same class
    for kW=1:nWrongPer
        theta1 = theta0;
        while (norm(theta1-theta0)<1e-8)
            k=ceil(rand*nTrial);
            theta1 = AllParameters(ParInds,k)';
        end
        ClassifData_XC((nWrongPer+1)*j-nWrongPer+kW,:) = [x0 theta1 0];
    end
end
ClassifData_XC(ClassifData_XC(:,1)==0,:)=[];

useNN=1; % Use regularization to prevent overfitting!
tic
[trainedClassifier, validationAccuracy] = ...
    trainClassifier(ClassifData_XC(:,1:end-1),...
    ClassifData_XC(:,end),useNN,nNodes,Lams(iL));
validationAccuracy
toc

save(strcat('TC',num2str(seed),'_XCorNoise',num2str(nlev(iL)),'PCA',num2str(SizesToTry(jEncSize)),'Lam',...
    num2str(Lams(iL)),'_nLay',num2str(nNodes),'_nW',num2str(nWrongPer),'ONEnP',num2str(nP),'.mat'), ...
   'trainedClassifier','numModes','validationAccuracy', '-v7.3');
end
end
end
end