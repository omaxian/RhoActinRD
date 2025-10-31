p = [0.55 0.35 0.3 4 1 1 0.05 1 0.1 10]; % Starfish
[rts,~,~,~] = PDERoots(p(1:9),0,20,100);
Stats=RhoAndActinPDEMod(p,0.1,1,1,3);
%X1=Stats;
    % InterpolatedSim=ResampleXCor(Stats.XCor,Stats.tSim,Stats.rSim,...
    %            ResampledX,ResampledT,10,120);
% Plan: try truncated Fourier representation (see how it works for different data 
% sets), and try data along a column (1 um spacing)