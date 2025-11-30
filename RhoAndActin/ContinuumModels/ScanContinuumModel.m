% koff0, rf, Nuc0, NucEn, koffAct, Arandom, Dv
% Generate bunch of samples in 2D
function []=ScanContinuumModel(RandomSeed)
addpath(genpath('..'))
try
    RandomSeed=str2num(RandomSeed);
catch
end
rng(RandomSeed);
nSamp = 250;
nP=10;
Params = zeros(nSamp,nP);
PBounds = [0.45 0.75; 0.01 0.5; 0 0.5; 0.5 7];
for iSamp=1:nSamp
    iSamp
    p=PBounds(:,1)+(PBounds(:,2)-PBounds(:,1)).*rand(size(PBounds,1),1);
    Nreg1D = round(20/p(4));
    Areg = 400/Nreg1D^2;
    Params(iSamp,:)=[p(1) 0.35 0.3*p(2) 4*p(2) p(2) p(3) 0.05 1 0.1 Areg];
    Stats=RhoAndActinDiffusionPDEs(Params(iSamp,:),0.1,1,1,...
        (RandomSeed-1)*nSamp+iSamp,0);
    AllStats{iSamp}=Stats;
    rng("shuffle")
end
save(strcat('Scan4P_',num2str(RandomSeed),'.mat'),'Params','AllStats');
end