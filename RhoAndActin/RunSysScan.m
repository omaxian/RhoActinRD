% For parameter scanning over the forward model
function []=RunSysScan(RandomSeed,nPar)
% Bounds for params
addpath(genpath('.'))
nSamp = 5;
nSeed = 10;
nParams = 8;
AllParameters=zeros(nParams,nSamp);
rng(RandomSeed);
for iSamp=1:nSamp
    iSamp
    if (iSamp>1)
        rng("shuffle")
    end
    xr=rand(5,1);
    if (nPar==5)
        Params = [0.7; 0.1+xr(1)*0.6; 1+xr(2)*39; 1; 1; ...
            0.1+xr(3)*9.9; 2; xr(4)*100];
        MonomerClock = Params(6)/Params(4)+Params(6)/Params(5)+Params(3);
        NucRate_B = Params(7)/(MonomerClock*Params(6));
        NucRate_Rho = Params(8)/(MonomerClock*Params(6));
        Params(7)=NucRate_B;
        Params(8)=NucRate_Rho;
    elseif (nPar==2)
        Params = [0.7; 0.4; 1+xr(2)*39; 1; 1; ...
            0.1+xr(3)*9.9; 0.054; 0.86];
    end
    AllParameters(:,iSamp)=Params;
    nNzs=zeros(nSeed,1);
    for seed=1:nSeed
        AllStats{seed}=RhoAndActinBasalNuc(Params,seed,0);
    end
    AllStatsAllSeeds{iSamp}=AllStats;
end
save(strcat('Scan',num2str(nPar),'_',num2str(RandomSeed),'.mat'),'AllStatsAllSeeds',...
    'AllParameters')
