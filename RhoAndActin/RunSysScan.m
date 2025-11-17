% For parameter scanning
function []=RunSysScan(RandomSeed)
% Bounds for params
addpath(genpath('.'))
nSamp = 200;
nSeed = 10;
nParams = 8;
numNonZero = 1;
AllParameters=zeros(nParams,nSamp);
AllnNzs = zeros(nSamp,1);
rng(str2num(RandomSeed));
for iSamp=1:nSamp
    iSamp
    if (iSamp>1)
        rng("shuffle")
    end
    xr=rand(5,1);
    Params = [0.5; 0.4; 1+xr(2)*39; 1; 1; ...
        0.1+xr(3)*9.9; xr(4)*7; xr(5)*50];
    MonomerClock = Params(6)/Params(4)+Params(6)/Params(5)+Params(3);
    NucRate_B = Params(7)/(MonomerClock*Params(6));
    NucRate_Rho = Params(8)/(MonomerClock*Params(6));
    Params(7)=NucRate_B;
    Params(8)=NucRate_Rho;
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
        end
        AvgStats.MeanRhoHat = AvgStats.MeanRhoHat/nNzTot;
        AvgStats.MeanActinHat = AvgStats.MeanActinHat/nNzTot;
        AvgStats.ACorsRho = AvgStats.ACorsRho/nNzTot;
        AvgStats.ACorsAct = AvgStats.ACorsAct/nNzTot;
        AvgStats.XCor = AvgStats.XCor/nNzTot;
        AvgStats.MeanActin = AvgStats.MeanActin/nNzTot;
    end
    AllAveragedStats{iSamp}=AvgStats;
end
save(strcat('ScanFullDataAuto_',num2str(RandomSeed),'.mat'),'AllAveragedStats',...
    'AllParameters','AllnNzs')
