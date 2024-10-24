% Compare the cross correlation generated from simulation to that from
% experiments
Bement=ZeroEr<20; % set to 1 for frog and starfish, 0 for C elegans
% Evaluate the cross correlation function
padxy=0;
% Cut out first 40 s
AllActin=AllActin(:,:,41:end);
AllRho=AllRho(:,:,41:end);
% Compute average size of Rho excitations
try
Thres=StSt(2);
catch
Thres=0.5;
end
RhoThres=AllRho>Thres;
[~,~,nFr]=size(RhoThres);
NumExcitations=zeros(nFr,1);
AvgExcitationSize=zeros(nFr,1);
for iT=1:nFr
    CC = bwconncomp(RhoThres(:,:,iT));
    L2=CC2periodic(CC,[1 1],'L');
    NumExcitations(iT)=max(L2(:));
    ThisExSize = zeros(NumExcitations(iT),1);
    for iJ=1:max(L2(:))
        ThisExSize(iJ)=sum(L2(:)==iJ)/(Nx^2)*L^2;
    end
    AvgExcitationSize(iT)=mean(ThisExSize);
end
[rSim,tSim,XCorsSim] = CrossCorrelations(dx,dx,dt*saveEvery,...
    AllRho,AllActin,padxy);
% Load the experimental
if (Bement)
    load('BementXCorsDS.mat')
    % Limit experimental to reasonable range
    tmax=120;
    rmax=10;
else
    load('BaixueXCorsDS.mat')
    % Shorter-ranged in C elegans
    tmax=60;
    rmax=5;
end
DistsByR=DistsByR(abs(dtvals)<=tmax, Uvals<=rmax);
Uvals=Uvals(Uvals<=rmax);
dtvals=dtvals(abs(dtvals)<tmax);
XCorsSim=XCorsSim(abs(tSim)<=tmax, rSim<=rmax);
rSim=rSim(rSim<=rmax);
tSim=tSim(abs(tSim)<=tmax);
% Interpolate the simulation to the experimental
InterpolatedSim = zeros(size(DistsByR));
for iT=1:length(dtvals)
    ttarg = dtvals(iT);
    % Find the window
    tdiff = abs(tSim-ttarg);
    [xs, index] = sort(tdiff);
    Inds=index(1:2);
    t1=tSim(Inds(1));
    t2=tSim(Inds(2));
    if ((t1-ttarg)*(t2-ttarg)>0)
        % Just take the closest one
        w=1;
    else
        w=(ttarg-t2)/(t1-t2); % is the weight for t1
    end
    if (w<0.5)
        keyboard
    end
    SpatialAtTime = w*XCorsSim(Inds(1),:)+(1-w)*XCorsSim(Inds(2),:);
    % Now reorganize in space
    for dR=1:length(Uvals)
        rtarg=Uvals(dR);
        % Find the window
        rdiff = abs(rSim-rtarg);
        [xs, index] = sort(rdiff);
        Inds=index(1:2);
        r1=rSim(Inds(1));
        r2=rSim(Inds(2));
        if ((r1-rtarg)*(r2-rtarg)>0)
            % Just take the closest one
            w=1;
        else
            w=(rtarg-r2)/(r1-r2); % is the weight for r1
        end
        if (w<0.5)
            keyboard
        end
        InterpolatedSim(iT,dR)=w*SpatialAtTime(Inds(1))...
            +(1-w)*SpatialAtTime(Inds(2));
    end
end
% Do norm of the difference
if (max(abs(InterpolatedSim(:)>0)))
    InterpolatedSim=InterpolatedSim/max(abs(InterpolatedSim(:)));
end
DistsByR=DistsByR/max(abs(DistsByR(:)));
% Distance weighted norms
DistsWts=exp(-Uvals'/2);
Difference=DistsWts.*(InterpolatedSim-DistsByR).*(InterpolatedSim-DistsByR);
% tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
% nexttile
% imagesc(rSim,tSim,XCorsSim/max(abs(XCorsSim(:))))
% clim([-1 1])
% nexttile
% imagesc(Uvals,dtvals,DistsByR)
% clim([-1 1])
% colormap turbo
% pause(1)
DiffNorm=sum(Difference(:))

