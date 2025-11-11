% For parameter scanning
function []=RunSysScan(RandomSeed)
% Bounds for params
addpath('Inputs/')
nSamp = 200;
nSeed = 10;
nParams = 8;
numNonZero = 2;
AllParameters=zeros(nParams,nSamp);
AllnNzs = zeros(nSamp,1);
rng(str2num(RandomSeed));
for iSamp=1:nSamp
    iSamp
    if (iSamp>1)
        rng("shuffle")
    end
    xr=rand(5,1);
    Params = [0.5+(xr(1)< 0.5)*0.2; 0.4; xr(2)*40; 1; 1; ...
        xr(3)*10; xr(4)*0.3; xr(5)*3];
    AllParameters(:,iSamp)=Params;
    nNzs=zeros(nSeed,1);
    for seed=1:nSeed
        AllStats{seed}=RhoAndActinBasalNuc(Params,seed,0);
        if (AllStats{seed}.XCor(1)~=0)
            nNzs(seed)=1;
        else
            nNzs(seed)=0;
        end
        if (sum(nNzs)==numNonZero)
            break;
        end
    end
    GoodInds = find(nNzs==1);
    nNzTot = length(GoodInds);
    AllnNzs(iSamp)=nNzTot;
    AvgStats = [];
    if (nNzTot>0)
        AllStats=AllStats(GoodInds);
        AvgStats=AllStats{1};
        for iP=2:nNzTot
            AvgStats.MeanRhoHat = AvgStats.MeanRhoHat+AllStats{iP}.MeanRhoHat;
            AvgStats.MeanActinHat = AvgStats.MeanActinHat+AllStats{iP}.MeanActinHat;
            AvgStats.ACorsRho = AvgStats.ACorsRho+AllStats{iP}.ACorsRho;
            AvgStats.ACorsAct = AvgStats.ACorsAct+AllStats{iP}.ACorsAct;
            AvgStats.XCor = AvgStats.XCor+AllStats{iP}.XCor;
            AvgStats.MeanActin = AvgStats.MeanActin+AllStats{iP}.MeanActin;
            AvgStats.ExSizes = [AvgStats.ExSizes; AllStats{iP}.ExSizes];
            AvgStats.NumStims =  AvgStats.NumStims+AllStats{iP}.NumStims;
        end
        AvgStats.MeanRhoHat = AvgStats.MeanRhoHat/nNzTot;
        AvgStats.MeanActinHat = AvgStats.MeanActinHat/nNzTot;
        AvgStats.ACorsRho = AvgStats.ACorsRho/nNzTot;
        AvgStats.ACorsAct = AvgStats.ACorsAct/nNzTot;
        AvgStats.XCor = AvgStats.XCor/nNzTot;
        AvgStats.MeanActin = AvgStats.MeanActin/nNzTot;
        AvgStats.NumStims = AvgStats.NumStims/nNzTot;
    end
    AllAveragedStats{iSamp}=AvgStats;
end
save(strcat('ScanFullData_',num2str(RandomSeed),'.mat'),'AllAveragedStats',...
    'AllParameters','AllnNzs')
