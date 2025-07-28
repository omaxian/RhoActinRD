% Make the bifurcation diagram and simulate the ODEs
koffs=0:0.001:2;
Pts=zeros(3,length(koffs));
for iK=1:length(koffs)
    kp0=0.05;
    kfb=1;
    KFB=0.1;
    kgapr=koffs(iK);
    p=0:0.001:100;
    OnRate = (kp0+kfb*p.^3./(KFB+p.^3));
    OffRate=(kgapr*p);
    Net=OnRate-OffRate;
    SgnChg=find((Net(1:end-1).*Net(2:end))<0);
    Pts(:,iK)=[p(SgnChg); zeros(3-length(SgnChg),1)];
    %plot(p,OnRate-OffRate)
    %hold on
    %set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1)
    %scatter(p(SgnChg),OnRate(SgnChg)-OffRate(SgnChg),'filled')
end
sw1=find(sum(Pts>0)>1,1,'first');
sw2=find(sum(Pts>0)>1,1,'last');
Branch1=[Pts(1,2:sw1-1) Pts(3,sw1:sw2)];
Branch2=Pts(2,sw1:sw2);
Branch3=Pts(1,sw1:end);
figure;
plot(koffs(2:sw2),Branch1,'-k')
hold on
plot(koffs(sw1:sw2),Branch2,':k')
plot(koffs(sw1:end),Branch3,'-k')
ylim([0 3])

figure;
% Solution of ODEs
pTime = zeros(100000,1);
kgapr = 1.6;
pTime(1)=1.25;
dt=1e-4;
for k=1:length(pTime)
OnRate = (kp0+kfb*pTime(k).^3./(KFB+pTime(k).^3));
OffRate=(kgapr*pTime(k));
pTime(k+1)=pTime(k)+dt*(OnRate-OffRate);
end
hold on
plot((0:length(pTime)-1)*dt,pTime)