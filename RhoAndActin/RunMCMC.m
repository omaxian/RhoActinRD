nEnsemble = 8;
nSamp = 250; % samples per ensemble
nParams = 6;
load('SeedParamsNoDiff.mat')
% k0s=0.9;
% rfs=0.4;
% FullLifetime=15;
% GrShrnk = 0.1; % rate in um/s
% MaxLength = 1; % in um
% N0PerMax = 0.2;
% Params = [k0s;rfs;FullLifetime;GrShrnk;MaxLength;N0PerMax];
AllDifferencesModelExp=zeros(nSamp,nEnsemble);
AllParameters=zeros(nParams*nEnsemble,nSamp);
CurrentParams=SeedParams(:,1:nEnsemble);
CurrentParams=CurrentParams(:);
AllMeanActins=zeros(nSamp,nEnsemble);
Accepted=zeros(nSamp,nEnsemble);
LastAccept=ones(nEnsemble,1);
for iSamp=1:nSamp
    for k=1:nEnsemble
        if (iSamp>1)
            % Proposal
            j=k;
            while (j==k)
                j = ceil(rand*nEnsemble);
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
%             k0Step = 0.05*randn;
%             rfStep = 0.05*randn;
%             LifetimeStep =2*randn;
%             GrShrnkStep = 0.05*randn;
%             MaxLenStep = 0.2*randn;
%             N0PerMaxStep=0.04*randn; % symmetric proposal
%             DeltaParams=max([k0Step;rfStep;LifetimeStep;...
%                 GrShrnkStep;MaxLenStep;N0PerMaxStep],-Params);
%             Params = OldParams+DeltaParams;
        else
            Params = CurrentParams((k-1)*nParams+1:k*nParams);
        end
        [DiffNorm,InterpolatedSim,MeanActin]=RhoAndActin(Params,iSamp);
        AllDifferencesModelExp(iSamp,k)=DiffNorm;
        AllParameters((k-1)*nParams+1:k*nParams,iSamp)=Params;
        AllXCors{iSamp,k}=InterpolatedSim;
        AllMeanActins(iSamp,k)=MeanActin;
        % Calculate probability and accept/reject
        if (iSamp > 1)
            kT=3;
            Likelihood=exp(-AllDifferencesModelExp(iSamp,k)/kT);
            OldLikelihood=exp(-AllDifferencesModelExp(LastAccept(k),k)/kT);
            vPrior = Prior(MeanActin,Params(1));
            OldvPrior = Prior(AllMeanActins(LastAccept(k),k),...
                AllParameters(nParams*(k-1)+1,LastAccept(k)));
            pAcc=z^(nParams-1)*(Likelihood*vPrior)/(OldLikelihood*OldvPrior);
            if (rand < pAcc) % accept
                LastAccept(k)=iSamp;
                Accepted(iSamp,k)=1;
                CurrentParams((k-1)*nParams+1:k*nParams)=Params;
            end
        end
    end
end

function p = Prior(MeanActin,koff)
    p = exp(-2*(MeanActin-1)^2.*(MeanActin>1)).*...
        exp(-500*(MeanActin-0.2)^2.*(MeanActin<0.2)).*(koff<1.25);
end
