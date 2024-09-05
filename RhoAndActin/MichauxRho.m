koffs=0:0.01:2;
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
plot(koffs(2:sw2),Branch1,'-k')
hold on
plot(koffs(sw1:sw2),Branch2,':k')
plot(koffs(sw1:end),Branch3,'-k')
