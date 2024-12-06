% Load the cross correlation function and excitation distribution
addpath('Inputs/')
EmType = "nmy"; % OPTIONS: NMY,Ani-NMU,Cyk1,Starfish
ActinOnly = 1;
LoadExisting = 1;
if (LoadExisting)
    if (ActinOnly)
        load(strcat(EmType,'MCMCRun.mat'))
    else
        load(strcat(EmType,'MCMCRun_All.mat'))
    end
    SampStart=iSamp+1;
else
if (EmType=="Starfish")
    load('SortedParametersOnlyActin.mat') % Params
    load('BementXCorsDS.mat')    % Cross corr fcn
else
    if (ActinOnly)
        load('SortedParametersCEOnlyActin.mat')
    else
        load('SortedParametersCE.mat')
    end
    load(strcat(EmType,"_Input.mat"));
    XCorsExp=XCorFilt;
end
WtsByR = exp(-Uvals'/2);
WtsByT = exp(-abs(dtvals)'/120);
TotWts=WtsByR.*WtsByT;
XCorNorm=TotWts.*XCorsExp.^2;
ZeroEr = round(sum(XCorNorm(:)),1);
if (ActinOnly)
    nWalker = 20;
    nSamp=750;
else
    nWalker=50;
    nSamp=500;
end
nSeed = 5; % averages per parameter set
nParams = 6;
PBounds = [0.4 1.22; 0.55 1.5; 0 30; 0 5; 0 10; 0 1];
ParamsStart=AllParametersSort(:,1:nWalker);
CurrentParams=ParamsStart(:);
AllDiffNorms=zeros(nSeed,nSamp,nWalker);
AllParameters=zeros(nParams*nWalker,nSamp);
AllMeanActins=zeros(nSeed,nSamp,nWalker);
AllExSizeErs = zeros(nSeed,nSamp,nWalker);
Accepted=zeros(nSamp,nWalker);
LastAccept=ones(nWalker,1);
SampStart=1;
end
for iSamp=SampStart:nSamp
    iSamp
    for k=1:nWalker
        if (iSamp>1)
        Params = inf*ones(nParams,1);
        while (~inRange(Params,PBounds))
            % Proposal
            j=k;
            while (j==k)
                j = ceil(rand*nWalker);
            end
            % Sample from g
            a=1.2;
            z=0;
            while (z<1/a || z>a)
                u=rand;
                z=(u*(a^(1/2)-a^(-1/2)+a^(-1/2))).^2;
            end
            Xj = CurrentParams((j-1)*nParams+1:j*nParams);
            Xk = CurrentParams((k-1)*nParams+1:k*nParams);
            Params = Xj + z*(Xk-Xj);
            Params = max(Params,0);
        end
        else
            Params = CurrentParams((k-1)*nParams+1:k*nParams);
        end
        for seed=1:nSeed
            tic
            Stats=RhoAndActin(Params,seed);
            % Compute the norm relative to the experiment and the
            % difference in the excitation size (for C. elegans only)
            ExSizeDiff=0;
            if (EmType~="Starfish")
                if (Stats.XCor(1)==0)
                    ExSizeDiff=1;
                else
                    WtsEx=(dsHist/2:dsHist:400);
                    xp=histcounts(Stats.ExSizes,0:dsHist:400);
                    xp=xp/(sum(xp)*dsHist);
                    ExSizeDiff = sum((xp-SizeHist).*(xp-SizeHist).*WtsEx)...
                        /sum(SizeHist.*SizeHist.*WtsEx); %L^2 norm
                end
            end
            % Cross correlation difference
            XCorEr = 1;
            if (Stats.XCor(1)~=0) 
                InterpolatedSim=ResampleXCor(Stats.XCor,Stats.tSim,Stats.rSim,...
                    Uvals,dtvals,max(Uvals)+1e-3,max(dtvals)+1e-3);
                XCorEr = TotWts.*(InterpolatedSim-XCorsExp).^2;
                XCorEr = sum(XCorEr(:))/ZeroEr;
            end
            toc
            ForgetIt = XCorEr > 1;
            AllDiffNorms(seed,iSamp,k)=XCorEr;
            AllMeanActins(seed,iSamp,k)=Stats.MeanActin;
            AllExSizeErs(seed,iSamp,k)=ExSizeDiff;
            if (ForgetIt) % Throw out really bad parameter sets
                AllDiffNorms(seed+1:end,iSamp,k)=AllDiffNorms(seed,iSamp,k);
                AllMeanActins(seed+1:end,iSamp,k)=AllMeanActins(seed,iSamp,k);
                AllExSizeErs(seed+1:end,iSamp,k)=AllExSizeErs(seed,iSamp,k);
                break;
            end
        end
        AllParameters((k-1)*nParams+1:k*nParams,iSamp)=Params;
        % Calculate relevant statistics
        MeanEr = mean(AllDiffNorms(:,iSamp,k));
        MeanActin = mean(AllMeanActins(:,iSamp,k));
        MeanExSizeEr = mean(AllExSizeErs(:,iSamp,k));
        if (iSamp > 1)
            MeanErPrev = mean(AllDiffNorms(:,LastAccept(k),k));
            MeanActinPrev = mean(AllMeanActins(:,LastAccept(k),k));
            MeanExSizeErPrev = mean(AllExSizeErs(:,LastAccept(k),k));
            Ell = Likelihood(MeanEr,MeanExSizeEr);
            PrevEll = Likelihood(MeanErPrev,MeanExSizeErPrev);
            Pr = Prior(MeanEr,MeanActin);
            PrevPr = Prior(MeanErPrev,MeanActinPrev);
            pAcc=z^(nParams-1)*(Ell*Pr)/(PrevEll*PrevPr);
            if (rand < pAcc) % accept
                LastAccept(k)=iSamp;
                Accepted(iSamp,k)=1;
                CurrentParams((k-1)*nParams+1:k*nParams)=Params;
            end
        end
    end
    if (mod(iSamp,5)==0)
        if (ActinOnly)
            save(strcat(EmType,'MCMCRun.mat'))
        else
            save(strcat(EmType,'MCMCRun_All.mat'))
        end
    end
end
