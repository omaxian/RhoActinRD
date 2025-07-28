% This file will output the movies in Fig. 2, and provides a starting point
% for running the code with a given parameter set 
close all;
ps=load('ParamsForFig2.mat');   
nSeed=10;
numNonZero=2;

% Load in experimental data
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

%figure(1)
%tiledlayout(1,8,'Padding', 'none', 'TileSpacing', 'compact')
ParamsToRun=ps.ans;
for iType = 1:4
    ps=ParamsToRun(iType,:);
    ExSizesAll=[];
    nNz=0;
    TotActin=0;
    for seed=1:nSeed
        tic
        Stats=RhoAndActinBasalNuc(ps,seed,nNz==1); % use nNz==1 for Fig 2 in paper
        toc
        % Compute the norm relative to the experiment and the
        % difference in the excitation size (for C. elegans only)
        if (Stats.XCor(1)==0)
        else
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
            if (nNz==numNonZero)
                break
            end
        end
    end
    % Compute errors 
    XCorAvg=XCorAvg/nNz;
    InterpolatedSim=ResampleXCor(XCorAvg,tSimulated,rSimulated,...
                        AllUvals{iType},Alldtvals{iType},...
                        max(AllUvals{iType})+1e-3,...
                        max(Alldtvals{iType})+1e-3);
    XCorNorm=AllXCorsExp{iType}.^2;
    ZeroEr = round(sum(XCorNorm(:)),1);
    XCorEr = (InterpolatedSim-AllXCorsExp{iType}).^2;
    XCorEr = sum(XCorEr(:))/ZeroEr
    if (iType>1)
        xp=histcounts(ExSizesAll,0:dsHist:400);
        xp=xp/(sum(xp)*AlldsHist{iType});
        Wts = (dsHist/2:dsHist:400);
        ExSizeDiff = sum((xp-AllSizeHist{iType})...
            .*(xp-AllSizeHist{iType}).*Wts)...
            /sum(AllSizeHist{iType}.*AllSizeHist{iType}.*Wts) %L^2 norm
    else
        ExSizeDiff=0;
    end
    
    % Cross correlations compared to data
    % figure(1)
    % nexttile
    % imagesc(rSimulated,tSimulated,XCorAvg)
    % xlim([0 max(Uvals)])
    % ylim([-120 120])
    % clim([-1 1])
    % xlabel('$\Delta r$ ($\mu$m)')
    % 
    % nexttile
    % imagesc(AllUvals{iType},Alldtvals{iType},AllXCorsExp{iType})
    % xlim([0 max(Uvals)])
    % ylim([-120 120])
    % clim([-1 1])
    % colormap turbo
    % xlabel('$\Delta r$ ($\mu$m)')
    % yticklabels('')
    
    % Excitation sizes 
    % figure(2)
    % if (iType>1)
    %     plot(dsHist/2:dsHist:400,xp)
    %     hold on
    %     set(gca,'ColorOrderIndex',iType-1)
    %     plot(dsHist/2:dsHist:400,AllSizeHist{iType},':')
    %     xlabel('Ex size ($\mu$m$^2$)')
    %     ylabel('pdf')
    %     xlim([0 200])
    % end
end

