% koff0, rf, Nuc0, NucEn, koffAct, Arandom, Dv
% Generate bunch of samples in 2D
function []=GenerateData(RandomSeed)
addpath(genpath('..'))
try
    RandomSeed=str2num(RandomSeed);
catch
end
rng(RandomSeed);
nSamp = 250;
nP=10;
Params = zeros(nSamp,nP);
PBounds = [0.01 0.5; 0 0.5; 0.5 10];
for iSamp=1:nSamp
    iSamp
    p=PBounds(:,1)+(PBounds(:,2)-PBounds(:,1)).*rand(3,1);
    Nreg1D = round(20/p(3));
    Areg = 400/Nreg1D^2;
    Params(iSamp,:)=[0.68 0.35 0.3*p(1) 4*p(1) 1*p(1) p(2) 0.05 1 0.1 Areg];
    Stats=RhoAndActinPDEMod(Params(iSamp,:),0.1,1,1,(RandomSeed-1)*nSamp+iSamp,0);
    AllStats{iSamp}=Stats;
    rng("shuffle")
end
save(strcat('Scan3P_',num2str(RandomSeed),'.mat'),'Params','AllStats');
end