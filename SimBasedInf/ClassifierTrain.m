% Main file to train a classifier from the data
load('ScanRGAInPrior.mat');
load('PrincipalComponentsXCor.mat')
nTrial = size(AllParameters,2);
ParInds=[2 3 6 10 11]; % parameters that vary
nlev = [5e-3]; % noise level
Lams = [1e-5]; % L2 regularization strength
for iL=1:length(Lams)
SizesToTry = [40]; % number of PCA modes
for jEncSize=1:length(SizesToTry)
% Load data
XCorScores = U'*AllXCors; % These are the scores
numModes=SizesToTry(jEncSize);
XCorScores = XCorScores(1:numModes,:);
Noise = nlev(iL)*SingVals(1:numModes).*randn(numModes,nTrial);
XCorScores = XCorScores+Noise;

nPVary=length(ParInds);
ClassifData_XC=zeros(2*nTrial,nPVary+1+size(XCorScores,1));
Weights=ones(2*nTrial,1);
for j=1:nTrial
    theta0 = AllParameters(ParInds,j)';
    x0 = XCorScores(:,j)';
    % Some other parameter not theta0 and not in same class
    k=j;
    while (k==j)
        k=ceil(rand*nTrial);
    end
    theta1 = AllParameters(ParInds,k)';
    ClassifData_XC(2*j-1:2*j,:) = [x0 theta0 1; x0 theta1 0];
end
Weights(2:2:end)=Weights(2:2:end)/mean(Weights(2:2:end));

useNN=1; % Use regularization to prevent overfitting!
tic
[trainedClassifier, validationAccuracy] = ...
    trainClassifier(ClassifData_XC(:,1:end-1),...
    ClassifData_XC(:,end),useNN,Weights,Lams(iL));
validationAccuracy
toc

save(strcat('TC_XCorNoise',num2str(nlev(iL)),'PCA',num2str(SizesToTry(jEncSize)),'Lam',...
    num2str(Lams(iL)),'_Hyb5.mat'), ...
   'trainedClassifier','numModes','validationAccuracy', '-v7.3');
end
end
