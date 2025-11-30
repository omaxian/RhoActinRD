function AllParams=AddParams(AllParams)
    MonomerClock = AllParams(:,6)./AllParams(:,4)+AllParams(:,3);
    AllParams(:,9)=MonomerClock;
    TotalLifetime = AllParams(:,6)./AllParams(:,4)+...
        AllParams(:,6)./AllParams(:,5)+AllParams(:,3);
    NFilBulk = AllParams(:,7).*TotalLifetime;
    MeanBulk = NFilBulk.*AllParams(:,6);
    AllParams(:,10)=MeanBulk;
    NFilRho = AllParams(:,8).*TotalLifetime;
    MeanRho = NFilRho.*AllParams(:,6);
    AllParams(:,11)=MeanRho;
end