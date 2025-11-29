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
PBounds = [0 0.5; 0 5; 0.01 0.5; 0 0.5; 0.5 7];
for iSamp=1:nSamp
    iSamp
    c=0;
    while (~c)
    p=PBounds(:,1)+(PBounds(:,2)-PBounds(:,1)).*rand(5,1);
    Nreg1D = round(20/p(5));
    Areg = 400/Nreg1D^2;
    Params(iSamp,:)=[0.68 0.35 p(1) p(2) p(3) p(4) 0.05 1 0.1 Areg];
    [rts,~,~,~] = PDERoots(Params(iSamp,1:9),0,0,0);
    c = length(rts(:,1))>1; % Has to diffuse more than 1/2 grid cell
    end
    Stats=RhoAndActinPDEMod(Params(iSamp,:),0.1,1,1,(RandomSeed-1)*nSamp+iSamp,0);
    AllStats{iSamp}=Stats;
    rng("shuffle")
end
save(strcat('Scan3P_',num2str(RandomSeed),'.mat'),'Params','AllStats');
end