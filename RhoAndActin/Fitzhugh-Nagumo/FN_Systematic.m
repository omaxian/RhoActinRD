% % Data
% load('BementXCorsDS.mat')    % Cross corr fcn
% XCorsExp{1}=DistsByR;
% tExp{1}=dtvals;
% rExp{1}=Uvals;
% XCorNorm=XCorsExp{1}.^2;
% ZeroErs(1) = round(sum(XCorNorm(:)),1);
% Names=["nmy" "nmy-pfn" "nmy-cyk"];
% for iName=2:4
%     load(strcat(Names(iName-1),"_Input.mat"));
%     XCorsExp{iName}=XCorFilt;
%     tExp{iName}=dtvals;
%     rExp{iName}=Uvals;
%     SizeHists{iName}=SizeHist;
%     XCorNorm=XCorsExp{iName}.^2;
%     ZeroErs(iName) = round(sum(XCorNorm(:)),1);
%     SizeNorms(iName) = sum(SizeHists{iName}.*SizeHists{iName});
% end
% % Get list of parameters
% FNSteadyStates;
% DActins=[0 1e-3 1e-2 0.1 1 10];
% nSamp = length(Params)*length(DActins);
% numNonZero = 1; % averages per parameter set
% nSeed = 1; % maximum # of attempts to get to 2
% nParams = 5;
% % a, b, Iext, tau, Dv
% AllDiffNorms=ones(nSamp,4);
% AllParameters=zeros(nParams,nSamp);
% AllMeanActins=zeros(nSamp,1);
% AllExSizeErs = ones(nSamp,4);
% AllExSizeErsT = ones(nSamp,4);
% Nnzs = zeros(nSamp,1);
% AllBehaviors=zeros(nSamp,1);
% iSamp=0;
for iP=2000:length(Params)
    if (mod(iP,100)==0)
        iP
        save('SystematicScanFN.mat')
    end
    for iD=1:length(DActins)
        iSamp=iSamp+1;
        AllBehaviors(iSamp)=Behaviors(iP);
        P=[Params(:,iP);DActins(iD)];
        TotActin=0;
        dt=0.5;
        [Stats,uns]=FNDynamics(P,1,dt,0);
        while (uns)
            dt=dt/5;
            [Stats,uns]=FNDynamics(P,1,dt,0);
            if (dt<1e-2)
                break;
            end
        end
        % Pot-process
        ExSizesAll=[];
        if (Stats.XCor(1)~=0)
            for iCond=1:4
                InterpolatedSim=ResampleXCor(Stats.XCor,Stats.tSim,Stats.rSim,...
                        rExp{iCond},tExp{iCond},max(rExp{iCond})+1e-3,...
                        max(tExp{iCond})+1e-3);
                XCorEr = (InterpolatedSim-XCorsExp{iCond}).^2;
                AllDiffNorms(iSamp,iCond) = sum(XCorEr(:))/ZeroErs(iCond);
            end
            for iCond=2:4
                % Excitation size errors
                xp=histcounts(Stats.ExSizes,0:dsHist:400);
                xp=xp/(sum(xp)*dsHist);
                XpEr=xp-SizeHists{iCond};
                AllExSizeErs(iSamp,iCond) = sum(XpEr.*XpEr)/SizeNorms(iName); %L^2 norm
                if (isnan(AllExSizeErs(iSamp,iCond)))
                    AllExSizeErs(iSamp,iCond)=1;
                end
            end
            for iCond=2:4
                % Excitation size errors
                xp=histcounts(Stats.ExSizesStThres,0:dsHist:400);
                xp=xp/(sum(xp)*dsHist);
                XpEr=xp-SizeHists{iCond};
                AllExSizeErsT(iSamp,iCond) = sum(XpEr.*XpEr)/SizeNorms(iName); %L^2 norm
                if (isnan(AllExSizeErsT(iSamp,iCond)))
                    AllExSizeErsT(iSamp,iCond)=1;
                end
            end
        end
        AllParameters(:,iSamp)=P;
    end
end
save('SystematicScanFN.mat')