% % Extract cross correlations and excitation sizes from filtered data 
% tmax=120;
% rmax=5;
% Name='spd';
% numMovies=10;
% for MovieNum=1:numMovies
%     AnalyzeCorrelations;
%     clearvars -except AllExes DistsByRF dtvalsF NumExes UvalsF tmax rmax Name MovieNum numMovies
%     save(strcat("Processed_",Name,"_",num2str(MovieNum),".mat"));
% end
 
% Averaging and overall distribution
% Make all xcors the same size
MinT=inf;
isSD = zeros(1,numMovies);
for MovieNum=1:numMovies
    load(strcat("Processed_",Name,"_",num2str(MovieNum),".mat"));
    load(strcat(Name,"Info_",num2str(MovieNum),".mat"));
    if (length(dtvalsF)<MinT && Info.IsSD)
        MinT=length(dtvalsF);
    end
    isSD(MovieNum)=Info.IsSD;
end
dsHist=4;
% Only process SD movies
IndsToUse = find(isSD);
for MovieNum=IndsToUse
    load(strcat("Processed_",Name,"_",num2str(MovieNum),".mat"));
    Diff=length(dtvalsF)-MinT;
    dtvalsF=dtvalsF(Diff/2+1:end-Diff/2);
    DistsByRF=DistsByRF(Diff/2+1:end-Diff/2,:);
    xp=histcounts(AllExes,0:dsHist:400);
    xps(MovieNum,:)=xp/(sum(xp)*dsHist);
    if (MovieNum==IndsToUse(1))
        XCorFilt=DistsByRF/sum(isSD);
        ZeroPos=find(dtvalsF==0);
    else
        XCorFilt=XCorFilt+DistsByRF/sum(isSD);
    end
end
xps(~isSD,:)=[];
SizeHist = mean(xps);
Uvals=UvalsF;
dtvals=dtvalsF;
% Fit size hist with log
% xmarks=dsHist/2:dsHist:400;
% p=polyfit(xmarks(SizeHist>0),log(SizeHist(SizeHist>0)),1);
% NormalizedFit=exp(p(1)*xmarks)/(sum(exp(p(1)*xmarks))*dsHist);
% plot(xmarks,SizeHist,xmarks,NormalizedFit)
% SizeHistFit = NormalizedFit;
