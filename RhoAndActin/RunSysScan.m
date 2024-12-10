% Bounds for params
addpath('Inputs/')
EmType = "nmy-cyk"; % nmy, nmy-cyk, nmy-pfn
load(strcat(EmType,"_Input.mat"));
XCorsExp=XCorFilt;
WtsByR = exp(-Uvals'/2);
WtsByT = exp(-abs(dtvals)'/120);
TotWts=WtsByR.*WtsByT;
XCorNorm=TotWts.*XCorsExp.^2;
ZeroEr = round(sum(XCorNorm(:)),1);
nSamp = 2500;
nSeed = 5;
nParams = 6;
PBounds = [0.4 1.22; 0.55 1.5; 0 30; 0 5; 0 10; 0 1];
AllDiffNorms=zeros(nSeed,nSamp);
AllParameters=zeros(nParams,nSamp);
AllExSizeErs = zeros(nSeed,nSamp);
rng(1);
for iSamp=1:nSamp
    % Make some random params in box
    r1=rand(6,1);
    params=PBounds(:,1)+r1.*(PBounds(:,2)-PBounds(:,1));
%     % Adjustment for actin inhibition
    % NewBoundsKinh = [max(PBounds(2,1)-params(1),0.2) PBounds(2,2)-params(1)];
    % params(2) = NewBoundsKinh(1)+r1(2)*(NewBoundsKinh(2)-NewBoundsKinh(1));
    params(1) = 0.8;
    params(2) = 0.4;
    close all;
    for seed=1:nSeed
        Stats=RhoAndActin(params,seed);
        % The criterion for moving on is a larger norm than zero PLUS a
        % local max (-1) in the cross correlation at (0,0)
        XCorEr = 1;
        if (Stats.XCor(1)~=0) 
            InterpolatedSim=ResampleXCor(Stats.XCor,Stats.tSim,Stats.rSim,...
                Uvals,dtvals,max(Uvals)+1e-3,max(dtvals)+1e-3);
            XCorEr = TotWts.*(InterpolatedSim-XCorsExp).^2;
            XCorEr = sum(XCorEr(:))/ZeroEr;
            WtsEx=(dsHist/2:dsHist:400);
            xp=histcounts(Stats.ExSizes,0:dsHist:400);
            xp=xp/(sum(xp)*dsHist);
            ExSizeDiff = sum((xp-SizeHist).*(xp-SizeHist).*WtsEx)...
                /sum(SizeHist.*SizeHist.*WtsEx); %L^2 norm
        else
            XCorEr=1;
            ExSizeDiff=1;
        end
        ForgetIt = XCorEr > 1;
        AllDiffNorms(seed,iSamp)=XCorEr;
        AllExSizeErs(seed,iSamp)=ExSizeDiff;
        if (ForgetIt) % Throw out really bad parameter sets
            AllDiffNorms(seed+1:end,iSamp)=AllDiffNorms(seed,iSamp);
            AllExSizeErs(seed+1:end,iSamp)=AllExSizeErs(seed,iSamp);
            break;
        end
    end
    AllParameters(:,iSamp)=params;
end
save('SystematicScan.mat')