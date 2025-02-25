load('SystematicScanFN.mat')
ExSizeErs=min(AllExSizeErs,AllExSizeErsT);
j=2;
[vals,ind]=sort(AllDiffNorms(:,j));
sp=AllParameters(:,ind);
BehSort=AllBehaviors(ind);
goodInds=find(vals<1);
% tiledlayout(3,2,'Padding', 'none', 'TileSpacing', 'compact')
% for iP=1:5
%     nexttile
%     histogram(sp(iP,goodInds))
% end
% return
p=AllParameters(:,ind(151));
p(end)=1;
[rts,Stab,spiral]=findFNRoots(p(1),p(2),p(3),p(4),0,0,0);
dt=0.5;
[Stats,unst]=FNDynamics(p,1,dt,1);
while(unst)
    dt=dt/5;
    [Stats,unst]=FNDynamics(p,1,dt,1);
end
nexttile
imagesc(Stats.rSim,Stats.tSim,Stats.XCor)
colormap turbo
clim([-1 1])
xlim([0 10])
ylim([-120 120])
xlabel('$\Delta t$ (s)')
ylabel('$\Delta r$ (s)')