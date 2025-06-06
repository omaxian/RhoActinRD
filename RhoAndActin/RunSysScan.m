function []=RunSysScan(RandomSeed)
% Bounds for params
addpath('Inputs/')
nSamp = 2500;
nSeed = 10;
nParams = 8;
numNonZero = 2;
AllDiffNorms=zeros(nSamp,4)*inf;
AllParameters=zeros(nParams,nSamp);
AllExSizeErs = zeros(nSamp,4)*inf;
AllMeanActins = zeros(nSamp,1);
nNzs = zeros(nSamp,1);
rng(str2num(RandomSeed));
EmTypes = ["Starfish" "nmy" "nmy-pfn" "nmy-cyk"];
for iType=1:4
    try
        load(strcat(EmTypes(iType),"_Input.mat"));
        AllXCorsExp{iType}=XCorFilt;
        Alldtvals{iType}=dtvals;
        AllUvals{iType}=Uvals;
        AllSizeHist{iType}=SizeHist;
        AlldsHist{iType}=dsHist;
    catch
        load('BementXCorsDS.mat')
        AllXCorsExp{iType}=DistsByR;
        Alldtvals{iType}=dtvals;
        AllUvals{iType}=Uvals;
    end
end
for iSamp=1:nSamp
    ExSizesAll=[];
    nNz=0;
    TotActin=0;
    if (iSamp>1)
        rng("shuffle")
    end
    xr=rand(5,1);
    Params = [0.8; 0.4; xr(1)*40; xr(2)*5; xr(2)*5; ...
        xr(3)*15; xr(4)*0.5; xr(5)*3];
    for seed=1:nSeed
        Stats=RhoAndActinBasalNuc(Params,seed,0);
        % Compute the norm relative to the experiment and the
        % difference in the excitation size (for C. elegans only)
        if (Stats.XCor(1)~=0)
            ExSizesAll=[ExSizesAll;Stats.ExSizes];
        end
        % Cross correlation difference
        XCorEr = 1;
        if (Stats.XCor(1)~=0) 
            nNz=nNz+1;
            if (nNz==1)
                XCorAvg=Stats.XCor;
            else
                XCorAvg=XCorAvg+Stats.XCor;
            end
            TotActin=TotActin+Stats.MeanActin;
            tSimulated=Stats.tSim;
            rSimulated=Stats.rSim;
        end
        if (nNz==numNonZero)
            break;
        end
    end
    % Compute errors 
    if (nNz>0)
        XCorAvg=XCorAvg/nNz;
        for iType=1:4
            InterpolatedSim=ResampleXCor(XCorAvg,tSimulated,rSimulated,...
                    AllUvals{iType},Alldtvals{iType},...
                    max(AllUvals{iType})+1e-3,...
                    max(Alldtvals{iType})+1e-3);
            XCorNorm=AllXCorsExp{iType}.^2;
            ZeroEr = round(sum(XCorNorm(:)),1);
            XCorEr = (InterpolatedSim-AllXCorsExp{iType}).^2;
            XCorEr = sum(XCorEr(:))/ZeroEr;
            AllDiffNorms(iSamp,iType)=XCorEr;
            if (iType > 1)
                xp=histcounts(ExSizesAll,0:dsHist:400);
                xp=xp/(sum(xp)*AlldsHist{iType});
                ExSizeDiff = sum((xp-AllSizeHist{iType})...
                    .*(xp-AllSizeHist{iType}))...
                    /sum(AllSizeHist{iType}.*AllSizeHist{iType}); %L^2 norm
                AllExSizeErs(iSamp,iType)=ExSizeDiff;
            end
        end
        AllMeanActins(iSamp)=TotActin/nNz;
        nNzs(iSamp)=nNz;
    end
    % Compute errors for each embryo
    AllParameters(:,iSamp)=Params;
end
save(strcat('Scank8_',num2str(RandomSeed),'.mat'))
