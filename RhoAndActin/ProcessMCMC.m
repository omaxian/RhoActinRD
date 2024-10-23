%load('GW25000.mat')
%load('FirstAttempt.mat')
nP=nParams;
AllDifferencesModelExp = mean(AllDiffNorms);
TwoMeanActins = mean(AllMeanActins);
MaxDiff = max(AllDiffNorms);
inds=1:length(AllDifferencesModelExp);
%AllParametersES=AllParameters;
%AllDifferencesModelExp=AllDifferencesModelExp(:);
%Accepted=Accepted(:);
%AllParameters=[];
%for p=1:nWalker
%    AllParameters = [AllParameters AllParametersES((p-1)*nP+1:p*nP,:)];
%end
%inds=find(AllDifferencesModelExp>0 & AllDifferencesModelExp<8 &Accepted);
Ap2Inds=find(abs(AllDifferencesModelExp-ZeroEr)<1e-10);
inds=setdiff(inds,Ap2Inds);
Ap3Inds=inds(abs(AllDiffNorms(10,inds)-AllDiffNorms(9,inds)) < 1e-10);
inds=setdiff(inds,Ap3Inds);
Ap4Inds=[];%inds(AllParameters(3,inds)>6);
inds=setdiff(inds,Ap4Inds);
%AllParameters(6,:)=TwoMeanActins;
subplot(1,3,1)
plot([0.55 0.55],[0 3],':k')
hold on
set(gca,'ColorOrderIndex',1)
scatter(AllParameters(1,Ap2Inds),AllParameters(2,Ap2Inds),10,'d','filled')
hold on
scatter(AllParameters(1,Ap3Inds),AllParameters(2,Ap3Inds),10,'^','filled')
scatter(AllParameters(1,inds),AllParameters(2,inds),50,AllDifferencesModelExp(inds),...
    's','filled')
ylim([0.2 1.1])
xlim([0.4 1.22])
colormap(flipud(turbo))
box on
xlabel('$k_\textrm{off}^{(0)}$')
ylabel('$k_\textrm{inh}$')
subplot(1,3,2)
scatter(AllParameters(3,Ap2Inds),AllParameters(5,Ap2Inds),10,'d','filled')
hold on
scatter(AllParameters(3,Ap3Inds),AllParameters(5,Ap3Inds),10,'^','filled')
scatter(AllParameters(3,inds),AllParameters(5,inds),...
    50,AllDifferencesModelExp(inds),'s','filled')
ylabel('$\ell_\textrm{max}$')
xlabel('$T_\textrm{fil}$')
box on
subplot(1,3,3)
scatter(AllParameters(4,Ap2Inds),AllParameters(6,Ap2Inds),10,'d','filled')
hold on
scatter(AllParameters(4,Ap3Inds),AllParameters(6,Ap3Inds),10,'^','filled')
scatter(AllParameters(4,inds),AllParameters(6,inds),...
    50,AllDifferencesModelExp(inds),'s','filled')
box on
ylabel('$\bar q_0$')
xlabel('$\nu_p$')