% Bounds for params
nSamp = 1000;
nSeed = 10;
nParams = 6;
ZeroEr = 14.4;
PBounds = [0.4 1.22; 0.55 1.5; 1 30; 0.1 5; 0.25 10; 0 0.6];
AllDifferencesModelExp=zeros(nSeed,nSamp);
AllParameters=zeros(nParams,nSamp);
AllMeanActins=zeros(nSeed,nSamp);
rng(1);
for iSamp=1:nSamp
    % Make some random params in box
    r1=rand(6,1);
    params = PBounds(:,1)+r1.*(PBounds(:,2)-PBounds(:,1));
    % Adjustment for actin inhibition
    %NewBoundsKinh = [max(PBounds(2,1)-params(1),0.2) PBounds(2,2)-params(1)];
    %params(2) = NewBoundsKinh(1)+r1(2)*(NewBoundsKinh(2)-NewBoundsKinh(1));
    params(2) = 0.4;
    DiffNorms=zeros(nSeed,1);
    close all;
    for seed=1:nSeed
        Stats=RhoAndActin(params,seed,ZeroEr);
        % The criterion for moving on is a larger norm than zero PLUS a
        % local max (-1) in the cross correlation at (0,0)
        if (Stats.DiffNorm==14.4)
            % It's zero
            ForgetIt=0;
        else
            ZeroZeroCor = Stats.InterpolatedSim(27,1);
            ForgetIt = abs(ZeroZeroCor+1) < 1e-2 && Stats.DiffNorm > ZeroEr;
            if (abs(ZeroZeroCor+1) < 1e-2 && Stats.DiffNorm < ZeroEr)
                seed
                warning('-1 cross corr and positive er')
            end
        end
        AllStats{seed}=Stats;
        if (ForgetIt) % Throw out really bad parameter sets
            DiffNorms(seed+1:end)=DiffNorms(seed);
            break;
        end
    end
    AllDifferencesModelExp(:,iSamp)=DiffNorms;
    AllParameters(:,iSamp)=params;
end
