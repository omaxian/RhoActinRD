nSamp=10000;
nSeed=2;
numNonZero=1;
AllnRrts=zeros(nSamp,1);
AllnUnst=zeros(nSamp,1);
Allps=zeros(9,nSamp);
AllDiffNorms=zeros(nSamp,1);
AllMeanActins=zeros(nSamp,1);
AllExSizeErs=zeros(nSamp,1);
Nnzs = zeros(nSamp,1);
SimSeeds = zeros(nSamp,nSeed);
EmType="nmy";
load(strcat(EmType,"_Input.mat"));
XCorsExp=XCorFilt;
TotWts=ones(length(dtvals),length(Uvals));
XCorNorm=TotWts.*XCorsExp.^2;
ZeroEr = round(sum(XCorNorm(:)),1);

for iSamp=1:nSamp
    rng("shuffle");
    % 1 - (0.3,0.6) % 2 - (0,0.3) % 3 - (0,1) % 4 - (0,10) 
    % 5 - (0,0.5) % 6 - (0,5) 
    ps = [0.3+rand*0.3; rand*0.3; rand; 10*rand; 0.5*rand; 5*rand; ...
        0.05; 1; 0.1];
    [rts,stability] = PDERoots(ps,0,20,100);
    % Check the stability (don't accept anything with one unif state)
    OneUnst = isscalar(rts(:,1)) && stability(1)==-1;
    TwoUnst = sum(stability==-1)>1;
    Bistable = sum(stability==-1)==1 && length(rts(:,1))==3;
    Allps(:,iSamp)=ps;
    AllnRrts(iSamp)=size(rts,1);
    AllnUnst(iSamp)=sum(stability==-1);
    % Run the simulation
    nNz=0;
    TotActin=0;
    ExSizesAll=[];
    for seed=1:nSeed
        %inseed = rand*100000; % for reproducibility
        %SimSeeds(iSamp,seed)=inseed;
        % Run (only the first time) to find the stability limit
        [Stats,st]=RhoAndActinTauPDEs(ps,seed,1);
        if (EmType~="Starfish")
            if (Stats.XCor(1)==0)
            else
                ExSizesAll=[ExSizesAll;Stats.ExSizes];
            end
        end
        % Compute the norm relative to the experiment and the
        % difference in the excitation size (for C. elegans only)
        % Cross correlation difference
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
    InterpolatedSim=ResampleXCor(XCorAvg,tSimulated,rSimulated,...
                Uvals,dtvals,max(Uvals)+1e-3,max(dtvals)+1e-3);
    XCorEr = TotWts.*(InterpolatedSim-XCorsExp).^2;
    XCorEr = sum(XCorEr(:))/ZeroEr;
    if (EmType~="Starfish")
    xp=histcounts(ExSizesAll,0:dsHist:400);
    WtsEx=ones(1,length(xp));
    xp=xp/(sum(xp)*dsHist);
    ExSizeDiff = sum((xp-SizeHist).*(xp-SizeHist).*WtsEx)...
            /sum(SizeHist.*SizeHist.*WtsEx); %L^2 norm
    if (isnan(ExSizeDiff))
        ExSizeDiff=1;
    end
    else
    ExSizeDiff=0;
    end
    MActin=TotActin/nNz;
    else
    XCorEr=2;
    ExSizeDiff=2;
    MActin=0;
    end
    Nnzs(iSamp)=nNz;
    AllDiffNorms(iSamp)=XCorEr;
    AllMeanActins(iSamp)=MActin;
    AllExSizeErs(iSamp)=ExSizeDiff;
end

% Find the best samples 
[vals,inds]=sort(AllDiffNorms(1:iSamp-1));
nRtsSrt=AllnRrts(inds);
StabSrt=AllnUnst(inds);


