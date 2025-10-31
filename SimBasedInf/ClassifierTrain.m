% koff0, rf, Nuc0, NucEn, koffAct, Dv
% Load data
load('DataForClassif_TwoParamsk45.mat')
% Identify parameter sets that are not doing anything
MAbove=zeros(length(PctAboveSaddle),1);
for jM=1:length(PctAboveSaddle)
MAbove(jM)=mean(PctAboveSaddle{jM});
end
SittingThereHi = MAbove>0.99;
SittingThereLow = MAbove < 0.01;
AllXCors(:,SittingThereHi | SittingThereLow)=0;
AllACorsAct(:,SittingThereHi | SittingThereLow)=0;
AllACorsRho(:,SittingThereHi | SittingThereLow)=0;
[nSetsOG,~]=size(AllParams);
nSets = 10000;
Inds=randperm(nSetsOG,nSets);
AllParams=AllParams(Inds,:);
%xInds = (31:2:91)+(0:10)'*121;
ResampledT = -120:2:120;
ResampledX = 0:0.5:10;
[X,T]=meshgrid(ResampledX,ResampledT);
X0Inds=find(X==0 & T>=-60 & T<=60);
T0Inds=find(T==0);
xInds=unique([T0Inds(1:2:end);X0Inds(1:2:end)]);
%xInds=1:length(AllXCors(:,1));
AllXCors=AllXCors(xInds,Inds);
nM=10;
gInds=(1:nM)+(0:10:(nM-1)*10)';
gInds=gInds(2:10);
%gInds=gInds(2:end);
nC=5;
cInds=(1:nC)+(0:10:(nC-1)*10)';
cInds = [cInds];
AllMeanActHat=AllMeanActHat(gInds,Inds)./AllMeanActHat(1,Inds);
AllMeanRhoHat=AllMeanRhoHat(gInds,Inds)./AllMeanRhoHat(1,Inds);
AllACorsAct=AllACorsAct(cInds,Inds);
AllACorsRho=AllACorsRho(cInds,Inds);
SittingThereHi = SittingThereHi(Inds);
SittingThereLow = SittingThereLow(Inds);
ParInds=[5 6];
nPVary=length(ParInds);
ClassifData_XC=zeros(2*nSets,nPVary+1+size(AllXCors,1));
ClassifData_MeanF=zeros(2*nSets,nPVary+1+2*size(AllMeanRhoHat,1));
ClassifData_ACorF=zeros(2*nSets,nPVary+1+2*size(AllACorsAct,1));
ClassifData_MeanAndACor=zeros(2*nSets,nPVary+1+2*size(AllMeanRhoHat,1)+2*size(AllACorsAct,1));
for j=1:nSets
    theta0 = [AllParams(j,ParInds)];
    x0 = AllXCors(:,j)';
    f0 = [AllMeanRhoHat(:,j)' AllMeanActHat(:,j)'];
    s0 = [AllACorsRho(:,j)' AllACorsAct(:,j)'];
    % Some other parameter not theta0 and not in same class
    k=j;
    normdP = norm(AllParams(j,ParInds)./AllParams(k,ParInds)-1);
    while (normdP<1e-2 || (SittingThereHi(k)==1 && SittingThereHi(j)==1) || ...
            (SittingThereLow(k)==1 && SittingThereLow(j)==1))
        k=ceil(rand*nSets);
        theta1 = [AllParams(k,ParInds)];
        normdP = norm(AllParams(j,ParInds)./AllParams(k,ParInds)-1);
    end
    ClassifData_XC(2*j-1:2*j,:) = [x0 theta0 1; x0 theta1 0];
    ClassifData_MeanF(2*j-1:2*j,:) = [f0 theta0 1; f0 theta1 0];
    ClassifData_ACorF(2*j-1:2*j,:) = [s0 theta0 1; s0 theta1 0];
    ClassifData_MeanAndACor(2*j-1:2*j,:) = [f0 s0 theta0 1; f0 s0 theta1 0];
end
useNN=1;
tic
[trainedClassifier, validationAccuracy] = ...
    trainClassifier(ClassifData_MeanF(:,1:end-1),ClassifData_MeanF(:,end),useNN);
validationAccuracy
toc
save(strcat('TrainedClassif_RelMeanFour1D',num2str(nM),'_2P'), ...
    'trainedClassifier','gInds','cInds', '-v7.3');

tic
[trainedClassifier, validationAccuracy] = ...
    trainClassifier(ClassifData_XC(:,1:end-1),ClassifData_XC(:,end),useNN);
validationAccuracy
toc
save('TrainedClassif_XCor_2P', ...
    'trainedClassifier','xInds', '-v7.3');



