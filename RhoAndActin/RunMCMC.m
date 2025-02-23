% Load the cross correlation function and excitation distribution
addpath('Inputs/')
EmType = "Starfish"; % nmy, nmy-cyk, nmy-pfn, star
ActinOnly = 1;
LoadExisting = 1;
if (LoadExisting)
    if (ActinOnly)
        load(strcat(EmType,'MCMCRunTwoEx_NoWt.mat'))
        nSamp=2000;
    else
        load(strcat(EmType,'MCMCRunTwoEx_NoWt_All.mat'))
        nSamp=1000;
    end
    SampStart=iSamp+1;
else
if (EmType=="Starfish")
    %load('SortedParametersOnlyActin.mat') % Params
    load('BementXCorsDS.mat')    % Cross corr fcn
    XCorsExp=DistsByR;
else
    if (ActinOnly)
        load('SortedParametersCEOnlyActin.mat')
    else
        load('SortedParametersCE.mat')
    end
    load(strcat(EmType,"_Input.mat"));
    XCorsExp=XCorFilt;
end
%WtsByR = exp(-Uvals'/2);
%WtsByT = exp(-abs(dtvals)'/120);
%TotWts=WtsByR.*WtsByT;
TotWts=ones(length(dtvals),length(Uvals));
XCorNorm=TotWts.*XCorsExp.^2;
ZeroEr = round(sum(XCorNorm(:)),1);
if (ActinOnly)
    nWalker = 20;
    nParamsV=4;
else
    nWalker = 40;
    nParamsV = 6;
end
nSamp = 2000;
numNonZero = 2; % averages per parameter set
nSeed = 10; % maximum # of attempts to get to 2
nParams = 6;
PBounds = [0.4 1.22; 0.55 1.5; 0 30; 0 5; 0 20; 0 1];
%ParamsStart=AllParametersSort(:,1:nWalker);
%CurrentParams=ParamsStart(:);
AllDiffNorms=zeros(nSamp,nWalker);
AllParameters=zeros(nParams*nWalker,nSamp);
AllMeanActins=zeros(nSamp,nWalker);
AllExSizeErs = zeros(nSamp,nWalker);
Nnzs = zeros(nSamp,nWalker);
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
        ExSizesAll=[];
        nNz=0;
        TotActin=0;
        for seed=1:nSeed
            Stats=RhoAndActin(Params,seed);
            % Compute the norm relative to the experiment and the
            % difference in the excitation size (for C. elegans only)
            if (EmType~="Starfish")
                if (Stats.XCor(1)==0)
                else
                    ExSizesAll=[ExSizesAll;Stats.ExSizes];
                end
            end
            % Cross correlation difference
            XCorEr = 1;
            if (Stats.XCor(1)~=0) 
                nNz=nNz+1;
                if (nNz==1)
                    XCorAvg=Stats.XCor;
                else
                    XCorAvg=XCorAvg+Stats.XCor;
                end
                TotActin=TotActin+Stats.MeanActin;
                tSimulated=Stats.tSim;
                rSimulated=Stats.rSim;
            end
            if (nNz==numNonZero)
                break;
            end
        end
        % Compute errors 
        if (nNz>0)
        XCorAvg=XCorAvg/nNz;
        InterpolatedSim=ResampleXCor(XCorAvg,tSimulated,rSimulated,...
                    Uvals,dtvals,max(Uvals)+1e-3,max(dtvals)+1e-3);
        XCorEr = TotWts.*(InterpolatedSim-XCorsExp).^2;
        XCorEr = sum(XCorEr(:))/ZeroEr;
        if (EmType~="Starfish")
        xp=histcounts(ExSizesAll,0:dsHist:400);
        WtsEx=ones(1,length(xp));
        xp=xp/(sum(xp)*dsHist);
        ExSizeDiff = sum((xp-SizeHist).*(xp-SizeHist).*WtsEx)...
                /sum(SizeHist.*SizeHist.*WtsEx); %L^2 norm
        else
        ExSizeDiff=0;
        end
        MActin=TotActin/nNz;
        else
        XCorEr=1;
        ExSizeDiff=1;
        MActin=0;
        end
        Nnzs(iSamp,k)=nNz;
        AllDiffNorms(iSamp,k)=XCorEr;
        AllMeanActins(iSamp,k)=MActin;
        AllExSizeErs(iSamp,k)=ExSizeDiff;
        AllParameters((k-1)*nParams+1:k*nParams,iSamp)=Params;
        % Calculate relevant statistics
        MeanEr = AllDiffNorms(iSamp,k);
        MeanActin = AllMeanActins(iSamp,k);
        MeanExSizeEr = AllExSizeErs(iSamp,k);
        if (iSamp > 1)
            MeanErPrev = AllDiffNorms(LastAccept(k),k);
            MeanActinPrev = AllMeanActins(LastAccept(k),k);
            MeanExSizeErPrev = AllExSizeErs(LastAccept(k),k);
            Ell = Likelihood(MeanEr,MeanExSizeEr);
            PrevEll = Likelihood(MeanErPrev,MeanExSizeErPrev);
            Pr = Prior(MeanEr,MeanActin);
            PrevPr = Prior(MeanErPrev,MeanActinPrev);
            pAcc=z^(nParamsV-1)*(Ell*Pr)/(PrevEll*PrevPr);
            if (rand < pAcc) % accept
                LastAccept(k)=iSamp;
                Accepted(iSamp,k)=1;
                CurrentParams((k-1)*nParams+1:k*nParams)=Params;
            end
        end
    end
    if (mod(iSamp,5)==0)
        if (ActinOnly)
            save(strcat(EmType,'MCMCRunTwoEx_NoWt.mat'))
        else
            save(strcat(EmType,'MCMCRunTwoEx_NoWt_All.mat'))
        end
    end
end

function a = inRange(Params,Range)
    % Check first one and then modify second
    x1 = Params(1) > Range(1,1) && Params(1) < Range(1,2);
    if (~x1)
        a=0;
        return;
    end
    Range(2,:)=Range(2,:)-Params(1);
    Range(2,1)=max(Range(2,1),0.2);
    x1 = Params > Range(:,1) & Params < Range(:,2);
    a = sum(x1)==length(Params);
end

function ell = Likelihood(Ebar,ExSizeEr)
    kT = 0.08;
    ell=exp(-(Ebar+ExSizeEr)/kT);
end

function p = Prior(MeanEr,MeanActin)
    % Penalize mean actin < 0.2 and high errors
    p = exp(-200*(MeanActin-0.2)^2*(MeanActin<0.2))...
       *exp(-50*(MeanEr-1)*(MeanEr>1));
end