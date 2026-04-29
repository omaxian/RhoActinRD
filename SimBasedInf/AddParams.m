% This function takes the 8 parameters that are actual 
% inputs to the simulation (koff0,kGAP,Tfil,nup,nud,Lmax,qb,qRho), 
% and computes the subunit turnover time (Monomer Clock), and the amount
% of basal and Rho-mediated F-actin, appending these three extra variables
% to the input array
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