% This file makes Fig. 6A or B. It finds the one parameter that 
% could vary with kGAP, then plots the marginal probability distributions.
% Parameters
nDose = 5;
ParInds=[2 3 6 10 11]; % all varying parameters
VariedParams = [1 2]; % Index in ParInds; 1 = kGAP, 2 = Tfil, 4 = qb
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

N=50; % Number for the 3 parameters to scan over
PBounds = [0.2 0.6; 1 60; 0.1 10; 0 8; 0 200];
if (0)
% Varied params
pv = zeros(N,length(VariedParams));
for k=1:length(VariedParams)
    pv(:,k) = PBounds(VariedParams(k),1)+(1/2:N)./N.*...
        (PBounds(VariedParams(k),2)-PBounds(VariedParams(k),1));
end
InPar = zeros(N^2,length(ParInds));
[p1,p2]=meshgrid(pv(:,1),pv(:,2));
InPar(:,VariedParams(1))=p1(:);
InPar(:,VariedParams(2))=p2(:);
GoodSets=[];
LocMax=0;
for iP1=0:N-1
    NewPBounds=PBounds;
    if (any(VariedParams==2))
        % Find the length 
        Len = PBounds(FixedParams(1),1)+iP1./(N-1).*...
                (PBounds(FixedParams(1),2)-PBounds(FixedParams(1),1));
        MaxLifetime = PBounds(2,2)+2*Len;
        MaxBasalNucRate = PBounds(end-1,2)/(Len*MaxLifetime);
        MaxRhoNucRate = PBounds(end,2)/(Len*MaxLifetime);
        NewPBounds(end-1,2)=MaxBasalNucRate;
        NewPBounds(end,2)=MaxRhoNucRate;
    elseif (any(VariedParams==3))
        % Find the time 
        Lifetime = PBounds(FixedParams(1),1)+iP1./(N-1).*...
                (PBounds(FixedParams(1),2)-PBounds(FixedParams(1),1));
        MaxLifetime = 2*PBounds(3,2)+Lifetime;
        MaxBasalNucRate = PBounds(end-1,2)/(PBounds(3,2)*MaxLifetime);
        MaxRhoNucRate = PBounds(end,2)/(PBounds(3,2)*MaxLifetime);
        NewPBounds(end-1,2)=MaxBasalNucRate;
        NewPBounds(end,2)=MaxRhoNucRate;
    end
    for iP2 = 0:N-1
        for iP3 = 0:N-1
            Posterior = zeros(N^2,nDose);
            % Fixed parameter values
            pf = NewPBounds(FixedParams,1)+[iP1;iP2;iP3]./(N-1).*...
                (NewPBounds(FixedParams,2)-NewPBounds(FixedParams,1));
            InPar(:,FixedParams) = pf'.*ones(N^2,length(FixedParams));
            if (any(VariedParams==2) || any(VariedParams==3))
                TotalLifetime = 2*InPar(:,3)+InPar(:,2);
                NFilBulk = InPar(:,4).*TotalLifetime;
                InPar(:,4) = NFilBulk.*InPar(:,3);
                NFilRho = InPar(:,5).*TotalLifetime;
                InPar(:,5) = NFilRho.*InPar(:,3);
            end
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
                AvgPostkGAP = zeros(N,nDose);
                for p=1:N
                    AvgPostkGAP(p,:)=mean(Posterior((p-1)*N+(1:N),:));
                end
                NormalizedAvg = AvgPostkGAP./sum(AvgPostkGAP);
                kGAPVals = InPar(1:N:end,1);
                FirstMoment = sum(kGAPVals.*NormalizedAvg);
                PenaltyTerm = (FirstMoment(2:end)-FirstMoment(1:end-1));
                % for iDose = 2:nDose
                %     % Penalize 
                %     PenaltyTerm(iDose-1) = -sum((NormalizedAvg(:,iDose) - NormalizedAvg(:,iDose-1)));%...
                %         %.*(2*(kGAPVals<FirstMoment(iDose-1))-1));
                % end
                if (min(PenaltyTerm)> 0)
                    GoodSets=[GoodSets;pf' mean(Posterior) -PenaltyTerm];
                end
            elseif (min(max(Posterior))>LocMax)
               LocMax=min(max(Posterior));
               pmax=pf;
            end
            
        end
    end
end
else
    load("GoodSets_Tfil.mat")
end
%%
TotalLoss = min(GoodSets(:,length(FixedParams)+(1:nDose)),[],2)-5*max(GoodSets(:,length(FixedParams)+nDose+(1:nDose-1)),[],2);
% Define the "best fixed params as those which minimize total loss"
[~,bestfixed]=max(TotalLoss);
%GoodSets=GoodSets(bestfixed,:);
%bestfixed=1;
N2D=50; % Number for the 2 parameters to scan over
% Varied params
pv = zeros(N2D,length(VariedParams));
for k=1:length(VariedParams)
    pv(:,k) = PBounds(VariedParams(k),1)+(1/2:N2D)./N2D.*...
        (PBounds(VariedParams(k),2)-PBounds(VariedParams(k),1));
end
nGood = size(GoodSets,1);
Posterior = zeros(N2D^2,nDose);
MaxPByDose = zeros(length(ParInds),nDose);
MaxValByDose = zeros(1,nDose);
for iG=1:nGood
pf = GoodSets(iG,1:length(FixedParams));
InPar = zeros(N2D^2,length(ParInds));
[p1,p2]=meshgrid(pv(:,1),pv(:,2));
InPar(:,VariedParams(1))=p1(:);
InPar(:,VariedParams(2))=p2(:);
InPar(:,FixedParams) = pf.*ones(N2D^2,length(FixedParams));
if (any(VariedParams==2) || any(VariedParams==3))
    TotalLifetime = 2*InPar(:,3)+InPar(:,2);
    NFilBulk = InPar(:,4).*TotalLifetime;
    InPar(:,4) = NFilBulk.*InPar(:,3);
    NFilRho = InPar(:,5).*TotalLifetime;
    InPar(:,5) = NFilRho.*InPar(:,3);
end
% Remove params not in prior
[~,prsc] = TC_Sustainable.predictFcn(InPar);
InPrior=prsc(:,2)>priorCutoff;
InPrior(InPar(:,4)>1.6./InPar(:,1))=0;
InPrior=logical(InPrior);
ParamsInPrior=InPar(InPrior,:);
% Fix the fixed params and extract the 2D marginals
AvgLik = zeros(size(InPar,1),nDose);
for iDose=1:nDose
    DataXCor=XCorsByDose(:,:,iDose);
    TestXCors=DataXCor(:);
    TestXCorsTr=U'*TestXCors;
    TestXCorsTr = TestXCorsTr(1:nModes,:);
    % Use classifier to evaluate likelihood to evidence
    InputVec = [TestXCorsTr'.*ones(size(ParamsInPrior,1),1)...
        ParamsInPrior];
    for seed=1:nSeed
        [~,scores] = ClassifBySeed{seed}.predictFcn(InputVec); 
        AvgLik(InPrior,iDose)=AvgLik(InPrior,iDose)+1/nSeed*(scores(:,2)./(1-scores(:,2)));
    end
end
Posterior(InPrior,:) = Posterior(InPrior,:)+1/nGood*AvgLik(InPrior,:);
if (iG==bestfixed)
    [MaxValByDose,indsmax]=max(AvgLik);
    MaxPByDose=InPar(indsmax,:);
end
end

TotalEvidence = sum(nGood*Posterior);

figure;
tiledlayout(1,nDose,'Padding', 'none', 'TileSpacing', 'compact');
Names = ["$k_\textrm{GAP}$" "$T_\textrm{fil}$" "$\ell_\textrm{max}$" ...
    "$f_b$" "$f_\rho$"];
for iDose=1:nDose
    nexttile(iDose)
    contourf(InPar(1:N2D:end,VariedParams(1)),InPar(1:N2D,VariedParams(2)),...
        reshape(log10(Posterior(:,iDose)),N2D,N2D),...
        -2:0.5:floor(max(reshape(log10(Posterior(:,iDose)),[],1))))
    hold on
%    scatter(pc(iDose,2),pc(iDose,10),40,'k','filled')
   % if (iDose==5)
    colorbar
    %end
    t=turbo(100);
    c=t(30:85,:);
    colormap(gca,c)
    xlabel(Names(VariedParams(1)))
    if (iDose==1)
    ylabel(Names(VariedParams(2)))
    end
    title(strcat(num2str(Doses(iDose)),' ng/\muL'),'interpreter','tex','Fontweight','normal')
    pbaspect([1 1 1])
    %ylim([0.08 4])
    %xlim([0.2 0.6])
    %clim([-2 1])
end
InPar=MaxPByDose;
nPFill=size(InPar,1);
pMaxAll = [0.7*ones(nPFill,1) InPar(:,1:2) 1*ones(nPFill,2) InPar(:,3:5)];
MonomerClock = pMaxAll(:,6)./pMaxAll(:,4)+pMaxAll(:,6)./pMaxAll(:,5)+pMaxAll(:,3);
NucRate_B = pMaxAll(:,7)./(MonomerClock.*pMaxAll(:,6));
NucRate_Rho = pMaxAll(:,8)./(MonomerClock.*pMaxAll(:,6));
pMaxAll(:,7)=NucRate_B;
pMaxAll(:,8)=NucRate_Rho;
pMaxAll=AddParams(pMaxAll);


%%
if (1)
load('PsPlot_Qb.mat')
N2=N2D;
AvgLik=zeros(N2^2,nDose);
for iDose=1:nDose
   % Generate 2D grid in fb, frho space
    fb = (1/2:N2)/N2*8;
    frho = (1/2:N2)/N2*200;
    [fb,frho]=meshgrid(fb,frho);
    ParamsInf = [pc(iDose,[2 3 6]).*ones(N^2,1) fb(:) frho(:)];
    [~,prsc] = TC_Sustainable.predictFcn(ParamsInf);
    InPrior=prsc(:,2)>priorCutoff;
    InPrior(ParamsInf(:,4)>1.6./ParamsInf(:,1))=0;
    InPrior=logical(InPrior);
    ParamsInPrior=ParamsInf(InPrior,:); 
    DataXCor=XCorsByDose(:,:,iDose);
    TestXCors=DataXCor(:);
    TestXCorsTr=U'*TestXCors;
    TestXCorsTr = TestXCorsTr(1:nModes,:);
    % Use classifier to evaluate likelihood to evidence
    InputVec = [TestXCorsTr'.*ones(size(ParamsInPrior,1),1)...
        ParamsInPrior];
    for seed=1:nSeed
        [~,scores] = ClassifBySeed{seed}.predictFcn(InputVec); 
        AvgLik(InPrior,iDose)=AvgLik(InPrior,iDose)+1/nSeed*(scores(:,2)./(1-scores(:,2)));
    end
end

figure;
tiledlayout(1,nDose,'Padding', 'none', 'TileSpacing', 'compact');
Names = [
    "$f_b$" "$f_\rho$"];
for iDose=1:nDose
    nexttile(iDose)
    contourf(ParamsInf(1:N2D:end,4),ParamsInf(1:N2D,5),...
        reshape(log10(AvgLik(:,iDose)),N2D,N2D),...
        -2:0.5:floor(max(reshape(log10(AvgLik(:,iDose)),[],1))))
    hold on
    load("PsPlot_Qb.mat")
    scatter(pc(iDose,10),pc(iDose,11),40,'k','filled')
    load("PsPlot_QbLowerQRho.mat")
    scatter(pc(iDose,10),pc(iDose,11),40,'b','filled')
    load("PsPlot_QbRaiseQRho.mat")
    purpleColor = [128, 0, 128] / 255; % Standard Purple RGB
    scatter(pc(iDose,10),pc(iDose,11),40,purpleColor,'filled')
   % if (iDose==5)
    colorbar
    %end
    t=turbo(100);
    c=t(30:85,:);
    colormap(gca,c)
    xlabel("$f_b$")
    if (iDose==1)
    ylabel("$f_\rho$")
    end
    title(strcat(num2str(Doses(iDose)),' ng/\muL'),'interpreter','tex','Fontweight','normal')
    pbaspect([1 1 1])
end
end