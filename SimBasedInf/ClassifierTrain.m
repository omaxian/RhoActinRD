Hybrid=0;
% Load data
if (Hybrid)
load('ScanUStim.mat')
else
load('Scan4PAll.mat') % continuum model
end
% Identify parameter sets that are not doing anything
if (Hybrid)
AllParameters=AddParams(AllParameters');
AllParameters = AllParameters';
end
nTrial = length(AllAveragedStats);
AllXCors=zeros(length(AllAveragedStats{1}.XCor(:)),nTrial);
AllMeanRhoHat=zeros(length(AllAveragedStats{1}.MeanRhoHat(:)),nTrial);
AllMeanActinHat=zeros(length(AllAveragedStats{1}.MeanActinHat(:)),nTrial);
AllACorsRho=zeros(length(AllAveragedStats{1}.ACorsRho(:)),nTrial);
AllACorsAct=zeros(length(AllAveragedStats{1}.ACorsAct(:)),nTrial);
AllInPrior=zeros(nTrial,1);
for j=1:nTrial
    AllXCors(:,j)=real(AllAveragedStats{j}.XCor(:));
    AllMeanRhoHat(:,j)=AllAveragedStats{j}.MeanRhoHat(:);
    AllMeanActinHat(:,j)=AllAveragedStats{j}.MeanActinHat(:);
    AllACorsRho(:,j)=AllAveragedStats{j}.ACorsRho(:);
    AllACorsAct(:,j)=AllAveragedStats{j}.ACorsAct(:);
    AllInPrior(j) = AllAveragedStats{j}.EnoughExcitation;
end

if (Hybrid)
ParInds=[3 6 10 11];
else
ParInds=[1 5 6 10];
end
% % Classifier for trajectories that sustain excitations (the prior)
useNN=1;
tic
[TC_Sustainable, validationAccuracy] = ...
    trainClassifier(AllParameters(ParInds,:)',AllInPrior>0,useNN);
toc
validationAccuracy
if (Hybrid)
save('TC_uStimInducSustain', 'TC_Sustainable', '-v7.3');
else
save('TC_ContSustain', 'TC_Sustainable', '-v7.3');
end


% Now filter the parameters to be only those that sustain excitations (as
% predicted by the classifier
InPrior = TC_Sustainable.predictFcn(AllParameters(ParInds,:)');
AllParameters=AllParameters(:,InPrior);
AllXCors=AllXCors(:,InPrior);
AllMeanRhoHat = AllMeanRhoHat(:,InPrior);
AllMeanActinHat = AllMeanActinHat(:,InPrior);
AllACorsRho=AllACorsRho(:,InPrior);
AllACorsAct=AllACorsAct(:,InPrior);
%xInds = (31:2:91)+(0:10)'*121;
%%
nTrial = size(AllParameters,2);
ResampledT = AllAveragedStats{1}.tSim;
ResampledX = AllAveragedStats{1}.rSim;
[X,T]=meshgrid(ResampledX,ResampledT);
X0Inds=find(X==0);% & mod(T,10)==0 & (T < 0 | T > 30));
T0Inds=find(T==0);%[];%find(T==-10 & X<=5 & mod(X,1)==0);
xInds=unique([T0Inds(1:2:end); X0Inds(1:2:end)]);
xInds=unique([T0Inds; X0Inds]);
if (Hybrid)
nM=5;
gInds=(1:nM)+(0:10:(nM-1)*10)';
gInds=gInds(2:end);
cInds=(201:205)';
else
nM=10;
gInds=(1:nM)+(0:10:(nM-1)*10)';
cInds=200+gInds(1:10);
gInds=gInds(2:10);
end
AllMeanActinHat=AllMeanActinHat(gInds,:)./AllMeanActinHat(1,:);
AllMeanRhoHat=AllMeanRhoHat(gInds,:)./AllMeanRhoHat(1,:);
% This is the autocorrelation at t=2
AllACorsAct=AllACorsAct(cInds,:);
AllACorsRho=AllACorsRho(cInds,:);
AllXCors=AllXCors(xInds,:);
nPVary=length(ParInds);
ClassifData_XC=zeros(2*nTrial,nPVary+1+size(AllXCors,1));
ClassifData_Four=zeros(2*nTrial,nPVary+1+2*size(AllMeanRhoHat,1)+2*size(AllACorsRho,1));
for j=1:nTrial
    theta0 = AllParameters(ParInds,j)';
    x0 = AllXCors(:,j)';
    f0 = [AllMeanRhoHat(:,j)' AllMeanActinHat(:,j)' AllACorsRho(:,j)' ...
        AllACorsAct(:,j)'];
    % Some other parameter not theta0 and not in same class
    k=j;
    while (k==j)
        k=ceil(rand*nTrial);
        theta1 = AllParameters(ParInds,k)';
    end
    ClassifData_XC(2*j-1:2*j,:) = [x0 theta0 1; x0 theta1 0];
    ClassifData_Four(2*j-1:2*j,:) = [f0 theta0 1; f0 theta1 0];
end

useNN=1;
tic
[trainedClassifier, validationAccuracy] = ...
    trainClassifier(ClassifData_XC(:,1:end-1),ClassifData_XC(:,end),useNN);
validationAccuracy
toc
if (Hybrid)
save('TC_XCor_UStimInducHyb', ...
   'trainedClassifier','xInds', '-v7.3');
else
save('TC_XCor_Cont', ...
   'trainedClassifier','xInds', '-v7.3');   
end

tic
[trainedClassifier, validationAccuracy] = ...
    trainClassifier(ClassifData_Four(:,1:end-1),ClassifData_Four(:,end),useNN);
validationAccuracy
toc
if (Hybrid)
save('TC_Four_UStimInducHyb', ...
   'trainedClassifier','gInds','cInds', '-v7.3');
else
save('TC_Four_Cont', ...
   'trainedClassifier','gInds','cInds', '-v7.3');
end



