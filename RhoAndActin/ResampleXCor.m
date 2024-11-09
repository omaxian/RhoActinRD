function InterpolatedSim=ResampleXCor(XCorsSim,tSim,rSim,Uvals,dtvals,rmax,tmax)
    Uvals=Uvals(Uvals<=rmax);
    dtvals=dtvals(abs(dtvals)<tmax);
    XCorsSim=XCorsSim(abs(tSim)<=tmax, rSim<=rmax);
    rSim=rSim(rSim<=rmax);
    tSim=tSim(abs(tSim)<=tmax);
    % Interpolate the simulation to the experimental
    InterpolatedSim = zeros(length(dtvals),length(Uvals));
    for iT=1:length(dtvals)
        ttarg = dtvals(iT);
        % Find the window
        tdiff = abs(tSim-ttarg);
        [~, index] = sort(tdiff);
        Inds=index(1:2);
        t1=tSim(Inds(1));
        t2=tSim(Inds(2));
        if ((t1-ttarg)*(t2-ttarg)>0)
            % Just take the closest one
            w=1;
        else
            w=(ttarg-t2)/(t1-t2); % is the weight for t1
        end
        SpatialAtTime = w*XCorsSim(Inds(1),:)+(1-w)*XCorsSim(Inds(2),:);
        % Now reorganize in space
        for dR=1:length(Uvals)
            rtarg=Uvals(dR);
            % Find the window
            rdiff = abs(rSim-rtarg);
            [~, index] = sort(rdiff);
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
end

