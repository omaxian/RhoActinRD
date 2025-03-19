% Load the cross correlation function and excitation distribution
addpath("../Inputs")
EmType = "nmy"; % nmy, nmy-cyk, nmy-pfn, star
LoadExisting = 0;
Randomness = 0;
if (LoadExisting)
    if (Randomness)
        load(strcat(EmType,'MCMCRunPDE_RndNuc.mat'))
    else
        load(strcat(EmType,'MCMCRunPDE_Det.mat'))
    end
    SampStart=iSamp+1;
    nSamp=1000;
else
if (EmType=="Starfish")
    %load('SortedParametersOnlyActin.mat') % Params
    load('BementXCorsDS.mat')    % Cross corr fcn
    XCorsExp=DistsByR;
else
    load(strcat(EmType,"_Input.mat"));
    XCorsExp=XCorFilt;
end
TotWts=ones(length(dtvals),length(Uvals));
XCorNorm=TotWts.*XCorsExp.^2;
ZeroEr = round(sum(XCorNorm(:)),1);
nSamp = 2000;
if (Randomness)
    nParams = 10;
    nParamsEff = 7; % number of actual varied params
    load('PartitionRegions.mat')
    nSeed = 2;
    PBounds = [0.3 0.6; 0 0.3; 0 1; 0 10; 0 0.5; 0 1; 0 1; ...
    0 inf; 0 inf; 0 inf];
else
    nParams = 9;
    nParamsEff = 6; % number of actual varied params
    nSeed = 1;
    PBounds = [0.3 0.6; 0 0.3; 0 1; 0 10; 0 0.5; 0 1; ...
    0 inf; 0 inf; 0 inf];
end
nWalker = 50;
%ParamsStart=AllParametersSort(:,1:nWalker);
%CurrentParams=P(:);
AllDiffNorms=zeros(nSamp,nWalker);
AllParameters=zeros(nParams*nWalker,nSamp);
AllMeanActins=zeros(nSamp,nWalker);
AllExSizeErs = zeros(nSamp,nWalker);
SimSeeds = zeros(nWalker*nSeed,nSamp);
Accepted=zeros(nSamp,nWalker);
LastAccept=ones(nWalker,1);
SampStart=1;
warning('Check bound on Df!')
end
for iSamp=SampStart:nSamp
    % Shuffle rng so you don't get same step every time
    iSamp
    for k=1:nWalker
        rng("shuffle");
        if (iSamp>1)
        ParamsRange=false;
        while (~ParamsRange)
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
            ParamsRange=inRange(Params,PBounds);
        end
        else
            Params = CurrentParams((k-1)*nParams+1:k*nParams);
        end
        ExSizesAll=[];
        TotActin=0;
        for seed=1:nSeed
            dt=0.1;
            if (Randomness)
                [Stats,st]=RhoAndActinPDEs_RandomNuc...
                    (Params,dt,seed,1,Nuc0s,NucEns);
            else
                [Stats,st]=RhoAndActinPDEs(Params,dt,1);
            end
            % Compute the norm relative to the experiment and the
            % difference in the excitation size (for C. elegans only)
            if (EmType~="Starfish")
                if (Stats.XCor(1)==0)
                else
                    ExSizesAll=[ExSizesAll;Stats.ExSizes];
                end
            end
            % Cross correlation difference
            if (Stats.XCor(1)~=0) 
                if (seed==1)
                    XCorAvg=Stats.XCor;
                else
                    XCorAvg=XCorAvg+Stats.XCor;
                end
                TotActin=TotActin+Stats.MeanActin;
                tSimulated=Stats.tSim;
                rSimulated=Stats.rSim;
            end
        end
        % Compute errors 
        XCorAvg=XCorAvg/nSeed;
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
            if (isnan(ExSizeDiff))
                ExSizeDiff=1;
            end
        else
            ExSizeDiff=0;
        end
        MActin=TotActin/nSeed;
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
            pAcc=z^(nParamsEff-1)*(Ell*Pr)/(PrevEll*PrevPr);
            if (rand < pAcc) % accept
                LastAccept(k)=iSamp;
                Accepted(iSamp,k)=1;
                CurrentParams((k-1)*nParams+1:k*nParams)=Params;
            end
        end
    end
    if (mod(iSamp,5)==0)
        if (Randomness)
            save(strcat(EmType,'MCMCRunPDE_RndNuc.mat'))
        else
            save(strcat(EmType,'MCMCRunPDE_Det.mat'))
        end
    end
end

function a = inRange(Params,Range)
    % Check first one and then modify second
    % x1 = Params(1) > Range(1,1) && Params(1) < Range(1,2);
    % if (~x1)
    %     a=0;
    %     return;
    % end
    % Range(2,:)=Range(2,:)-Params(1);
    % Range(2,1)=max(Range(2,1),0.2);
    x1 = Params > Range(:,1) & Params < Range(:,2);
    a = sum(x1)==length(Params);
    % Check the stability (don't accept anything with one unif state)
    % if (a==1)
    %     [rts,stability] = PDERoots(Params,0.1,20,100);
    %     %OneUnst = isscalar(rts(:,1)) && stability(1)==-1;
    %     %TwoUnst = sum(stability==-1)>1;
    %     OneSt = isscalar(rts(:,1)) && stability==1;
    %     if (OneSt)
    %         a=0;
    %     end
    % end
end

function ell = Likelihood(Ebar,ExSizeEr)
    kT = 0.08;
    ell=exp(-(Ebar+ExSizeEr)/kT);
end

function p = Prior(MeanEr,MeanActin)
    % Penalize mean actin < 0.2 and high errors
    p = 1;%exp(-500*(MeanActin-0.1)^2*(MeanActin<0.1));
end