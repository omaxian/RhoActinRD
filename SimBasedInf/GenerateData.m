% koff0, rf, Nuc0, NucEn, koffAct, Arandom, Dv
% Generate bunch of samples in 2D
function []=GenerateData(RandomSeed)
addpath(genpath('/home/ondrej/RhoModel/RhoActinRD'))
try
    RandomSeed=str2num(RandomSeed);
catch
end
rng(RandomSeed);
nSamp = 100;
nP=10;
Params = zeros(nSamp,nP);
% PBounds = [0.3 0.6; 0 0.3; 0 1; 0 15; 0 0.5; 0 1;...
%     0 inf; 0 inf; 0 inf; 1 50];
PBounds = [0.01 0.5; 0 0.5; 0.5 10];
for iSamp=1:nSamp
    %rts=[0 0];
    %while (isscalar(rts(:,1)))
    p=PBounds(:,1)+(PBounds(:,2)-PBounds(:,1)).*rand(3,1);
    Nreg1D = round(20/p(3));
    %Areg = 400/Nreg1D^2;
    Areg=(10/3)^2;
    Params(iSamp,:)=[0.45 0.45 0.3*p(1) 4*p(1) 1*p(1) p(2) 0.05 1 0.1 Areg];
    % Check for bistable possibilities
    %[rts,~,~,~] = PDERoots(Params(iSamp,1:9),1,20,100);
    %end
    Stats=RhoAndActinPDEMod(Params(iSamp,:),0.1,1,1,(RandomSeed-1)*nSamp+iSamp,0);
    AllStats{iSamp}=Stats;
    rng("shuffle")
end
save(strcat('Scan2PNew_',num2str(RandomSeed),'.mat'));
end