nSamp=1000000;
nSeed=5;
numNonZero=2;
AllnRrts=zeros(nSamp,1);
AllnUnst=zeros(nSamp,1);
Allps=zeros(9,nSamp);
AllDiffNorms=zeros(nSamp,1);
AllMeanActins=zeros(nSamp,1);
Nnzs = zeros(nSamp,1);
SimSeeds = zeros(nSamp,nSeed);
load('BementXCorsDS.mat')    % Cross corr fcn
XCorsExp=DistsByR;
TotWts=ones(length(dtvals),length(Uvals));
XCorNorm=TotWts.*XCorsExp.^2;
ZeroEr = round(sum(XCorNorm(:)),1);

for iSamp=1:nSamp
    ChangeMe=1;
    while (ChangeMe)
        % 1 - (0,1) % 2 - (0,5) % 3 - (0,1) % 4 - (Nuc0,5) % 5 - (0.4,1.5)
        % 6 - (0,0.2) % 7 - (0,0.6) % 8 - (2,20  % 9 - (0,10)
        ps = [rand*2; rand*5; rand; 10*rand; 2*rand; 0.1*rand; ...
            0.05; 1; 0.1];
        [rts,stability] = PDERoots(ps,0.1,20,100);
        % Check the stability (don't accept anything with one unif state)
        OneUnst = isscalar(rts(:,1)) && stability(1)==-1;
        TwoUnst = sum(stability==-1)>1;
        Bistable = sum(stability==-1)==1 && length(rts(:,1))==3;
        if (TwoUnst || OneUnst)
            ChangeMe=0;
        end
    end
    Allps(:,iSamp)=ps;
    AllnRrts(iSamp)=size(rts,1);
    AllnUnst(iSamp)=sum(stability==-1);
    % Run the simulation
    nNz=0;
    TotActin=0;
    for seed=1:nSeed
        %inseed = rand*100000; % for reproducibility
        %SimSeeds(iSamp,seed)=inseed;
        % Run (only the first time) to find the stability limit
        if (seed==1)
        st=0;
        dt=0.2;
        while (st==0)
            [~,st]=RhoAndActinPDEs(ps,dt,[],0);
            dt=dt/2;
            if (dt<0.001)
                break;
            end
        end
        end
        [Stats,st]=RhoAndActinPDEs(ps,dt,[],1);
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
    MActin=TotActin/nNz;
    else
    XCorEr=2;
    ExSizeDiff=2;
    MActin=0;
    end
    Nnzs(iSamp)=nNz;
    AllDiffNorms(iSamp)=XCorEr;
    AllMeanActins(iSamp)=MActin;
end

% Find the best samples 
[vals,inds]=sort(AllDiffNorms(1:iSamp-1));
nRtsSrt=AllnRrts(inds);
StabSrt=AllnUnst(inds);


