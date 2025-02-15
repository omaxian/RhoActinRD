nSamp=1000000;
AllnRrts=zeros(nSamp,1);
AllSt=zeros(nSamp,3);
AllSp=zeros(nSamp,3);
Allps=zeros(8,nSamp);
for j=1:nSamp
    ps=2*rand(8,1);
    Allps(:,j)=ps;
    [rts,stability,spiral] = PDERoots(ps);
    AllnRrts(j)=size(rts,1);
    AllSt(j,:)=stability;
    AllSp(j,:)=spiral;
end