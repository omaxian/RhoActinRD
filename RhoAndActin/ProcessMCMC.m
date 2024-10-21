%AllDifferencesModelExp=DifferencesModelExp;
%load('GoodmanWeareNoDiffBox.mat')
%load('GW25000.mat')
load('FirstAttempt.mat')
nP=nParams;
AllDifferencesModelExp = mean(AllDifferencesModelExp);
inds=1:length(AllDifferencesModelExp);
%AllParametersES=AllParameters;
%AllDifferencesModelExp=AllDifferencesModelExp(:);
%Accepted=Accepted(:);
%AllParameters=[];
%for p=1:nWalker
%    AllParameters = [AllParameters AllParametersES((p-1)*nP+1:p*nP,:)];
%end
%inds=find(AllDifferencesModelExp>0 & AllDifferencesModelExp<8 &Accepted);
Ap2Inds=find(abs(AllDifferencesModelExp-14.4)<1e-10);
inds=setdiff(inds,Ap2Inds);
Ap3Inds=find(AllDifferencesModelExp>14.4+1e-10);%inds(AllParameters(4,inds)>3.5 & AllParameters(3,inds)<5);
inds=setdiff(inds,Ap3Inds);
Ap4Inds=[];%inds(AllParameters(3,inds)>6);
inds=setdiff(inds,Ap4Inds);
subplot(1,3,1)
% plot([0.55 0.55],[0 3],':k')
% hold on
% set(gca,'ColorOrderIndex',1)
%scatter(AllParameters(1,:),AllParameters(2,:),10)
%hold on
% set(gca,'ColorOrderIndex',1)
% scatter(AllParameters(1,Accepted==1),AllParameters(2,Accepted==1)',10,'filled')
scatter(AllParameters(1,Ap2Inds),AllParameters(2,Ap2Inds),10,'d','filled')
hold on
scatter(AllParameters(1,Ap3Inds),AllParameters(2,Ap3Inds),20,'^','filled')
%scatter(AllParameters(1,Ap4Inds),AllParameters(2,Ap4Inds),50,'>','filled')
scatter(AllParameters(1,inds),AllParameters(2,inds),50,AllDifferencesModelExp(inds),...
    's','filled')
colormap(flipud(summer))
box on
xlabel('$k_\textrm{off}^{(0)}$')
ylabel('$k_\textrm{inh}$')
subplot(1,3,2)
%scatter(AllParameters(4,:),AllParameters(3,:),10)
%hold on
% set(gca,'ColorOrderIndex',1)
% scatter(AllParameters(4,Accepted==1),AllParameters(3,Accepted==1)',10,'filled')
scatter(AllParameters(4,Ap2Inds),AllParameters(3,Ap2Inds),10,'d','filled')
hold on
scatter(AllParameters(4,Ap3Inds),AllParameters(3,Ap3Inds),20,'^','filled')
scatter(AllParameters(4,inds),AllParameters(3,inds),50,AllDifferencesModelExp(inds),...
    's','filled')
%scatter(AllParameters(4,Ap4Inds),AllParameters(3,Ap4Inds),50,'>','filled')
xlabel('$\nu_p$')
ylabel('$T_\textrm{fil}$')
box on
subplot(1,3,3)
%scatter(AllParameters(6,:),AllParameters(5,:),10)
%hold on
% set(gca,'ColorOrderIndex',1)
% scatter(AllParameters(6,Accepted==1),AllParameters(5,Accepted==1)',10,'filled')
scatter(AllParameters(6,Ap2Inds),AllParameters(5,Ap2Inds),10,'d','filled')
hold on
scatter(AllParameters(6,Ap3Inds),AllParameters(5,Ap3Inds),20,'^','filled')
scatter(AllParameters(6,inds),AllParameters(5,inds),50,AllDifferencesModelExp(inds),...
    's','filled')
% scatter(AllParameters(6,Ap4Inds),AllParameters(5,Ap4Inds),50,'>','filled')
xlabel('$\bar q_\textrm{nuc}$')
ylabel('$\ell_\textrm{max}$')
box on
% subplot(1,4,4)
% scatter(AllParameters(7,Accepted==0),AllParameters(5,Accepted==0)')
% hold on
% set(gca,'ColorOrderIndex',1)
% scatter(AllParameters(7,Accepted==1),AllParameters(3,Accepted==1)','filled')
% scatter(AllParameters(7,inds),AllParameters(3,inds),50,'s','filled')
% scatter(AllParameters(7,Ap2Inds),AllParameters(3,Ap2Inds),50,'d','filled')
% scatter(AllParameters(7,Ap3Inds),AllParameters(3,Ap3Inds),50,'^','filled')
% scatter(AllParameters(7,Ap4Inds),AllParameters(4,Ap4Inds),50,'>','filled')
% xlabel('$D$')
% ylabel('$T_\textrm{fil}$')
% box on
