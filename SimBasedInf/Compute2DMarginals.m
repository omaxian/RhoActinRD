% Computes and plots the 2D likelihood function by averaging over
% other parameters. Inputs: likelihoodToevidence (array of values), 
% parametersInPrior, indices of varied parameters (ParInds), whether to
% make a plot, ParamsTest (if the true values are known), and some indices
% for the output figure (jD and nD).
% Outputs: the unique 1D values of all parameters, the 2D marginals as a
% set of matrices, and the 1D marginals over kGAP (if kGAP is varying).  
function [uvals,AllMarg,kraMarg] = Compute2DMarginals(LikelihoodToEvidence,ParametersInPrior,...
    ParInds,doPlot,ParamsTest,jD,nD)
    nV=length(ParInds);
    uvals = cell(nV,1);
    allindices = cell(nV,1);
    TotalNum=1;
    for k=1:nV
        [uvals{k},allindices{k}] = uniquetol(ParametersInPrior(:,ParInds(k)),1e-10,...
            'OutputAllIndices',true);
        TotalNum=TotalNum*length(uvals{k});
    end
    % Compute marginal
    if (length(ParInds)==4)
        IndsX = [2 3];
        IndsY = [1 4];
        xLabels = ["$\ell$" "$f_b$"];
        yLabels = ["$T_\textrm{fil}$" "$f_\rho$"];
        kraMarg=[];
    elseif (length(ParInds)==5)
        IndsX = [3 4];
        IndsY = [2 5];
        xLabels = ["$\ell$" "$f_b$"];
        yLabels = ["$T_\textrm{fil}$" "$f_\rho$"];
        % 1D marginal for kra
        if (length(uvals{1})>1)
        kraMarg = zeros(length(uvals{1}),1);
        nzation = TotalNum/length(uvals{1});
        for p=1:length(uvals{1})
            kraMarg(p) = sum(LikelihoodToEvidence(allindices{1}{p}))/nzation;
        end
        nexttile(jD)
        plot(uvals{1},kraMarg)
        xlabel('$k_\textrm{GAP}$')
        pbaspect([1 1 1])
        if (jD==1)
        ylabel('Mean LER','interpreter','tex')
        end
        xlim([0.2 0.6])
        end
        set(gca,'YScale','Log')
    end
    for jM=1:length(IndsX)
    ind1=IndsX(jM);
    ind2=IndsY(jM);
    FirstMarg = zeros(length(uvals{ind2}),length(uvals{ind1}));
    nzation = TotalNum/(length(uvals{ind1})*length(uvals{ind2}));
    for j=1:length(uvals{ind1})
        for k=1:length(uvals{ind2})
            inds = intersect(allindices{ind1}{j},allindices{ind2}{k});
            FirstMarg(k,j) = sum(LikelihoodToEvidence(inds))/nzation;
        end
    end
    AllMarg{jM}=FirstMarg;
    if (doPlot)
    nexttile%(jM*nD+jD)
    try
    contourf(uvals{ind1},uvals{ind2},log10(FirstMarg),-2:0.5:min(10,floor(max(log10(FirstMarg(:))))))
    catch
    %contourf(uvals{ind1},uvals{ind2},log10(FirstMarg))
    end
    %imagesc(uvals{ind1},uvals{ind2},log10(FirstMarg))
    set(gca,'YDir','Normal')
    colormap(gca,turbo)
    t=turbo(100);
    c=t(30:85,:);
    colormap(gca,c)
    %clim([-2 2])
    %colorbar
    try
        hold on
        scatter(ParamsTest(ParInds(ind1)),ParamsTest(ParInds(ind2)),50,'k','filled')
    catch
    end
    if (jD==nD && jM==length(IndsX))
        colorbar
    end
    if (jD==1)
        ylabel(yLabels(jM))
    end
    xlabel(xLabels(jM))
    pbaspect([1 1 1])
    end
    % figure(2)
    % Colors=get(gca,'ColorOrder');
    % GoodOnes = FirstMarg > 9;
    % % BW is your binary image
    % boundaries = bwboundaries(GoodOnes);
    % for k = 1:length(boundaries)
    %     boundary = boundaries{k};
    %     boundary(:,2) = uvals{ind1}(boundary(:,2));
    %     boundary(:,1) = uvals{ind2}(boundary(:,1));
    %     % Swap columns: boundary(:,2) is X, boundary(:,1) is Y
    %     fill(boundary(:,2), boundary(:,1),Colors(jD,:),'FaceAlpha', 0.25); 
    %     hold on
    %     plot(boundary(:,2), boundary(:,1),'Color',Colors(jD,:))
    %     hold on
    end
    % end
end

    

