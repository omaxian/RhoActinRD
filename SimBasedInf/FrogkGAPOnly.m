% This file is for Fig. S9 (considering the pssibility that only kGAP can
% vary). It is similar to OnePFrog.m
%function FrogkGAPOnly(iP1)
addpath(genpath('..'))
nDose = 5;
ParInds=[2 3 6 10 11]; % all varying parameters
VariedParams = [1]; % Index in ParInds; 1 = kGAP, 2 = Tfil, 4 = qb
FixedParams = setdiff(1:length(ParInds),VariedParams);
PostThres = 1; % LER threshold for a parameter set to be "good"
load('PrincipalComponentsXCor.mat')
load('XCorsByDose_Frog.mat')
% Load classifiers
nModes = 10;
priorCutoff=0.5;
Rzation=1e-3;
nSeed = 2;
nLay = [20 40 20];
nW=1;
load(strcat('TCPriorOne_5.mat'))
for seed=1:nSeed
    load(strcat('NewTC',num2str(seed-1),'_XCorNoise0',...
    'PCA',num2str(nModes),...
    'Lam',num2str(Rzation),'_nLay',num2str(nLay),...
    '_nW',num2str(nW),'ONEnP5.mat'),...
    'trainedClassifier')
    ClassifBySeed{seed}=trainedClassifier;
end

N=50; % Number of parameters to scan over
PBounds = [0.2 0.6; 1 60; 0.1 10; 0 8; 0 200];
% Varied params
pv = zeros(N,length(VariedParams));
for k=1:length(VariedParams)
    pv(:,k) = PBounds(VariedParams(k),1)+(1/2:N)./N.*...
        (PBounds(VariedParams(k),2)-PBounds(VariedParams(k),1));
end
InPar = pv;
GoodSets=[];
LocMax=0;
%iP1=iP1-1;
for iP1=0:N-1
    NewPBounds=PBounds;
    for iP2 = 0:N-1
        iP2
        for iP3 = 0:N-1
            for iP4 = 0:N-1
            Posterior = zeros(N,nDose);
            % Fixed parameter values
            pf = NewPBounds(FixedParams,1)+[iP1;iP2;iP3;iP4]./(N-1).*...
                (NewPBounds(FixedParams,2)-NewPBounds(FixedParams,1));
            pf = gs(bestfixed,1:4)';
            InPar(:,FixedParams) = pf'.*ones(N,length(FixedParams));
            % Remove params not in prior
            [~,priorsc] = TC_Sustainable.predictFcn(InPar);
            InPrior = priorsc(:,2) > priorCutoff;
            InPrior(InPar(:,4)>1.6./InPar(:,1))=0;
            InPrior=logical(InPrior);
            ParamsInPrior=InPar(InPrior,:);
            for iDose = 1:nDose
                rng(0);
                DataXCor=XCorsByDose(:,:,iDose);
                TestXCors=DataXCor(:);
                TestXCorsTr=U'*TestXCors;
                TestXCorsTr = TestXCorsTr(1:nModes,:);
                % Use classifier to evaluate likelihood to evidence
                InputVec = [TestXCorsTr'.*ones(size(ParamsInPrior,1),1)...
                    ParamsInPrior];
                for seed=1:nSeed
                    [~,scores] = ClassifBySeed{seed}.predictFcn(InputVec); 
                    Posterior(InPrior,iDose) = Posterior(InPrior,iDose)+1/nSeed*(scores(:,2)./(1-scores(:,2)));
                end
            end
            % Only keep parameter sets where max of posterior is good
            if (min(max(Posterior))>PostThres)
                % Compute average posterior for each kGAP
                NormalizedAvg = Posterior./sum(Posterior);
                kGAPVals = InPar(:,1);
                FirstMoment = sum(kGAPVals.*NormalizedAvg);
                PenaltyTerm = (FirstMoment(2:end)-FirstMoment(1:end-1));
                if (min(PenaltyTerm)> 0)
                    GoodSets=[GoodSets;pf' mean(Posterior) FirstMoment];
                end
            elseif (min(max(Posterior))>LocMax)
               LocMax=min(max(Posterior));
               pmax=pf;
            end
            end
        end
    end
%    save(strcat('Sets',num2str(iP1)),'GoodSets');
end