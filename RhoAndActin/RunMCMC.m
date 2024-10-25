Bement=1;
nWalker = 50;
nSamp = 500; % samples per ensemble
nSeed = 5; % averages per parameter set
nParams = 6;
PBounds = [0.4 1.22; 0.55 1.5; 0 30; 0 5; 0 10; 0 1];
if (Bement)
    ZeroEr=14.4;
    load('SortedParameters.mat')
    ParamsStart=AllParametersSort(:,1:nWalker);
    CurrentParams=ParamsStart(:);
else
    ZeroEr=261.6;
    load('SortedParametersCE.mat')
    ParamsStart=AllParametersSort(:,1:nWalker);
    CurrentParams=ParamsStart(:);
end
AllDiffNorms=zeros(nSeed,nSamp,nWalker);
AllParameters=zeros(nParams*nWalker,nSamp);
AllMeanActins=zeros(nSeed,nSamp,nWalker);
AllAverageExSize = zeros(nSeed,nSamp,nWalker);
AllNumEx = zeros(nSeed,nSamp,nWalker);
Accepted=zeros(nSamp,nWalker);
LastAccept=ones(nWalker,1);
for iSamp=1:nSamp
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
            Stats=RhoAndActin(Params,seed,ZeroEr);
            toc
            % The criterion for moving on is a larger norm than zero PLUS a
            % local max (-1) in the cross correlation at (0,0)
            if (abs(Stats.DiffNorm-ZeroEr) < 1e-5)
                ForgetIt=0; % It's zero
            else
                ForgetIt = Stats.DiffNorm > ZeroEr;
            end
            AllDiffNorms(seed,iSamp,k)=Stats.DiffNorm;
            AllMeanActins(seed,iSamp,k)=Stats.MeanActin;
            AllAverageExSize(seed,iSamp,k)=mean(Stats.AvgExcitation);
            AllNumEx(seed,iSamp,k)=mean(Stats.NumExcitations);
            if (ForgetIt) % Throw out really bad parameter sets
                AllDiffNorms(seed+1:end,iSamp,k)=AllDiffNorms(seed,iSamp,k);
                AllMeanActins(seed+1:end,iSamp,k)=AllMeanActins(seed,iSamp,k);
                AllAverageExSize(seed+1:end,iSamp,k)=AllAverageExSize(seed,iSamp,k);
                AllNumEx(seed+1:end,iSamp,k)=AllNumEx(seed,iSamp,k);
                break;
            end
        end
        AllParameters((k-1)*nParams+1:k*nParams,iSamp)=Params;
        % Calculate relevant statistics
        MeanEr = mean(AllDiffNorms(:,iSamp,k));
        MeanActin = mean(AllMeanActins(:,iSamp,k));
        MeanExSize = mean(AllAverageExSize(:,iSamp,k));
        if (iSamp > 1)
            MeanErPrev = mean(AllDiffNorms(:,LastAccept(k),k));
            MeanActinPrev = mean(AllMeanActins(:,LastAccept(k),k));
            MeanExSizePrev = mean(AllAverageExSize(:,LastAccept(k),k));
            Ell = Likelihood(MeanEr,ZeroEr);
            PrevEll = Likelihood(MeanErPrev,ZeroEr);
            Pr = Prior(MeanEr,MeanActin,MeanExSize,ZeroEr);
            PrevPr = Prior(MeanErPrev,MeanActinPrev,...
                MeanExSizePrev,ZeroEr);
            pAcc=z^(nParams-1)*(Ell*Pr)/(PrevEll*PrevPr);
            if (rand < pAcc) % accept
                LastAccept(k)=iSamp;
                Accepted(iSamp,k)=1;
                CurrentParams((k-1)*nParams+1:k*nParams)=Params;
            end
        end
    end
    if (mod(iSamp,10)==0)
        if (Bement)
            save('MCMCBementRun.mat')
        else
            save('MCMCCERun.mat')
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

function ell = Likelihood(Ebar,E0)
    kT=E0/7;
    ell=exp(-Ebar/kT);
end

function p = Prior(MeanEr,MeanActin,MeanExSize,ZeroEr)
    % Penalize mean actin < 0.2 and high errors
    p = exp(-200*(MeanActin-0.2)^2*(MeanActin<0.2))...
       *exp(-50*(MeanEr/ZeroEr-1)*(MeanEr>ZeroEr));
    if (ZeroEr > 20) % hack for Baixue
        % Add term for mean excitation size
        p=p*exp(-2e-4*(MeanExSize-20)^2);
    end 
end
