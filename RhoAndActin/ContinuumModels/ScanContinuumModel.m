% koff0, rf, Nuc0, NucEn, koffAct, Arandom, Dv
% Generate bunch of samples in 2D
function []=ScanContinuumModel(RandomSeed)
addpath(genpath('..'))
try
    RandomSeed=str2num(RandomSeed);
catch
end
rng(RandomSeed);
nSamp = 1000;
nP=10;
Params = zeros(nSamp,nP);
%PBounds = [0.45 0.75; 0.01 0.5; 0 0.5; 0.5 7];
%PBounds = [0 0.1; 0 0.5; 0.5 7];
%load('StarfishXCor.mat')
%load('WormXCor.mat')
for iSamp=1:nSamp
    % First choose the off rate for actin
    kdiss = 0.01+0.49*rand;
    Df = 0.5*rand;
    % Nucleation rates are proportional to diss rate
    qnuc0 = rand*2*kdiss;
    qnucRho = rand*20*kdiss;
    Lreg = 0.5+rand*(7-0.5);
    iSamp
    %p=PBounds(:,1)+(PBounds(:,2)-PBounds(:,1)).*rand(size(PBounds,1),1);
    Nreg1D = round(20/Lreg);
    Areg = 400/Nreg1D^2;
    %Params(iSamp,:)=[p(1) 0.35 0.3*p(2) 4*p(2) p(2) p(3) 0.05 1 0.1 Areg];
    Params(iSamp,:)=[0.35 0.35 qnuc0 qnucRho kdiss Df 0.05 1 0.1 Areg];
    Stats=RhoAndActinDiffusionPDEs(Params(iSamp,:),0.1,1,1,...
       (RandomSeed-1)*nSamp+iSamp,0);
    %ster=norm(DataXCor(:,1)-Stats.XCor(:,1))
   % Allsters(iSamp)=ster;
    AllStats{iSamp}=Stats;
    rng("shuffle")
end
save(strcat('Scan5P_',num2str(RandomSeed),'.mat'),'Params','AllStats');
end