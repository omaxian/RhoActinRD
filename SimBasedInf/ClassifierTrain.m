load('ScanUStimInPrior.mat');
load('PrincipalComponentsXCor.mat')
nTrial = size(AllParameters,2);
ParInds=[3 6 10 11];
nlev = [0 0 0 0 0 0];
Lams = [0 3.5e-5 1e-4 3.5e-4 1e-3 3.5e-3];% zeros(1,length(nlev));
for iL=1:length(Lams)
SizesToTry = [40];
for jEncSize=1:length(SizesToTry)
% Load data
% nTrial = size(AllParameters,2);
% DInds = randperm(nTrial,nData);
% AllParameters=AllParameters(:,DInds);
% AllXCors=AllXCors(:,DInds);

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
    if (0)
        x1=XCorScores(:,k)';
        RelDiff = abs(x1./x0-1)*SingVals(1:numModes);
        Weights(2*j)=RelDiff;
    end
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
    num2str(Lams(iL)),'_Hyb.mat'), ...
   'trainedClassifier','numModes','validationAccuracy', '-v7.3');
end
end
