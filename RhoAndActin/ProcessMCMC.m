%AllDifferencesModelExp=DifferencesModelExp;
load('GoodmanWeare.mat')
AllParametersES=AllParameters;
AllDifferencesModelExp=AllDifferencesModelExp(:);
Accepted=Accepted(:);
AllParameters=[];
for p=1:nEnsemble
    AllParameters = [AllParameters AllParametersES((p-1)*7+1:p*7,:)];
end
inds=find(AllDifferencesModelExp>0 & AllDifferencesModelExp<14 &Accepted);
Ap2Inds=inds(AllParameters(1,inds)<0.55 & AllParameters(2,inds)<1);
inds=setdiff(inds,Ap2Inds);
Ap3Inds=inds(AllParameters(3,inds) < 10 & AllParameters(4,inds)<2);
inds=setdiff(inds,Ap3Inds);
Ap4Inds=inds(AllParameters(1,inds)<0.55 & AllParameters(2,inds)>1);
inds=setdiff(inds,Ap4Inds);
subplot(1,4,1)
plot([0.55 0.55],[0 1.6],':k')
hold on
set(gca,'ColorOrderIndex',1)
scatter(AllParameters(1,Accepted==0),AllParameters(2,Accepted==0)')
set(gca,'ColorOrderIndex',1)
scatter(AllParameters(1,Accepted==1),AllParameters(2,Accepted==1)','filled')
scatter(AllParameters(1,inds),AllParameters(2,inds),50,'s','filled')
scatter(AllParameters(1,Ap2Inds),AllParameters(2,Ap2Inds),50,'d','filled')
scatter(AllParameters(1,Ap3Inds),AllParameters(2,Ap3Inds),50,'^','filled')
scatter(AllParameters(1,Ap4Inds),AllParameters(2,Ap4Inds),50,'>','filled')
xlabel('$k_\textrm{off}^{(0)}$')
ylabel('$k_\textrm{inh}$')
subplot(1,4,2)
scatter(AllParameters(4,Accepted==0),AllParameters(3,Accepted==0)')
hold on
set(gca,'ColorOrderIndex',1)
scatter(AllParameters(4,Accepted==1),AllParameters(3,Accepted==1)','filled')
scatter(AllParameters(4,inds),AllParameters(3,inds),50,'s','filled')
scatter(AllParameters(4,Ap2Inds),AllParameters(3,Ap2Inds),50,'d','filled')
scatter(AllParameters(4,Ap3Inds),AllParameters(3,Ap3Inds),50,'^','filled')
scatter(AllParameters(4,Ap4Inds),AllParameters(3,Ap4Inds),50,'>','filled')
xlabel('$\nu_p$')
ylabel('$T_\textrm{fil}$')
box on
subplot(1,4,3)
scatter(AllParameters(6,Accepted==0),AllParameters(5,Accepted==0)')
hold on
set(gca,'ColorOrderIndex',1)
scatter(AllParameters(6,Accepted==1),AllParameters(5,Accepted==1)','filled')
scatter(AllParameters(6,inds),AllParameters(5,inds),50,'s','filled')
scatter(AllParameters(6,Ap2Inds),AllParameters(5,Ap2Inds),50,'d','filled')
scatter(AllParameters(6,Ap3Inds),AllParameters(5,Ap3Inds),50,'^','filled')
scatter(AllParameters(6,Ap4Inds),AllParameters(5,Ap4Inds),50,'>','filled')
xlabel('$\bar q_\textrm{nuc}$')
ylabel('$\ell_\textrm{max}$')
box on
subplot(1,4,4)
scatter(AllParameters(7,Accepted==0),AllParameters(5,Accepted==0)')
hold on
set(gca,'ColorOrderIndex',1)
scatter(AllParameters(7,Accepted==1),AllParameters(3,Accepted==1)','filled')
scatter(AllParameters(7,inds),AllParameters(3,inds),50,'s','filled')
scatter(AllParameters(7,Ap2Inds),AllParameters(3,Ap2Inds),50,'d','filled')
scatter(AllParameters(7,Ap3Inds),AllParameters(3,Ap3Inds),50,'^','filled')
scatter(AllParameters(7,Ap4Inds),AllParameters(4,Ap4Inds),50,'>','filled')
xlabel('$D$')
ylabel('$T_\textrm{fil}$')
box on
