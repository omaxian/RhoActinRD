function [uvals,AllMarg] = Compute2DMarginals(LikelihoodToEvidence,ParametersInPrior,...
    ParInds,doPlot,ParamsTest,FixedParSet,jD,nD)
    if (FixedParSet==1)
        ParInds=ParInds(1:2);
    elseif (FixedParSet==2)
        ParInds=ParInds(end-1:end);
    end
    nV=length(ParInds);
    uvals = cell(nV,1);
    allindices = cell(nV,1);
    for k=1:nV
        [uvals{k},allindices{k}] = uniquetol(ParametersInPrior(:,ParInds(k)),1e-10,...
            'OutputAllIndices',true);
    end
    % Compute marginal
    if (FixedParSet==2)
        IndsX = [1];
        IndsY = [2];
        xLabels = ["$f_b$"];
        yLabels = ["$f_\rho$"];
    elseif (FixedParSet==1)
        IndsX = [2];
        IndsY = [1];
        xLabels = ["$\ell$"];
        yLabels = ["$T_\textrm{fil}$"];
    elseif (FixedParSet==0)
        IndsX = [2 3];
        IndsY = [1 4];
        xLabels = ["$\ell$" "$f_b$"];
        yLabels = ["$T_\textrm{fil}$" "$f_\rho$"];
    end
    for jM=1:length(IndsX)
    ind1=IndsX(jM);
    ind2=IndsY(jM);
    FirstMarg = zeros(length(uvals{ind2}),length(uvals{ind1}));
    for j=1:length(uvals{ind1})
        for k=1:length(uvals{ind2})
            inds = intersect(allindices{ind1}{j},allindices{ind2}{k});
            FirstMarg(k,j) = mean(LikelihoodToEvidence(inds));
        end
    end
    FirstMarg(isnan(FirstMarg))=0;
    AllMarg{jM}=FirstMarg;
    if (doPlot)
    nexttile%((jM-1)*nD+jD)
    try
    contourf(uvals{ind1},uvals{ind2},log10(FirstMarg),-5:0.5:min(10,floor(max(log10(FirstMarg(:))))))
    catch
    %contourf(uvals{ind1},uvals{ind2},log10(FirstMarg))
    end
    %imagesc(uvals{ind1},uvals{ind2},log10(FirstMarg))
    set(gca,'YDir','Normal')
    colormap(gca,turbo)
    %clim([-3 3])
    colorbar
    try
        hold on
        scatter(ParamsTest(ParInds(ind1)),ParamsTest(ParInds(ind2)),50,'k','filled')
    catch
    end
    ylabel(yLabels(jM))
    xlabel(xLabels(jM))
    pbaspect([1 1 1])
    end
    end
end

    

