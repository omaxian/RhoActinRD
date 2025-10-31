function IndChosen = ScatterPlotParams(LikelihoodToEvidence,...
    AllParams,ParamsTest,XCorsTest,ResampledX,ResampledT)
    [~,nData]=size(LikelihoodToEvidence);
    tiledlayout(3,3,...
        'Padding', 'none', 'TileSpacing', 'compact')
    IndChosen = zeros(nData,1);
    for jChk=1:nData
        [vals,inds]=sort(LikelihoodToEvidence(:,jChk),'ascend');
        valsp=vals;
        indsp=inds;
        nexttile%(2*(jChk-1)+1)
        % scatter3(AllParams(indsp,5),AllParams(indsp,6),sqrt(AllParams(indsp,10)),10,log10(valsp),'filled')
        % hold on
        % scatter3(ParamsTest(jChk,5),ParamsTest(jChk,6),sqrt(ParamsTest(jChk,10)),100,'ks','filled')
        % scatter3(AllParams(inds(end),5),AllParams(inds(end),6),sqrt(AllParams(inds(end),10)),100,'m^','filled')
        % %view(2)
        % ylabel('$k_\textrm{diss}$')
        % colormap jet
        % climlim=clim;
        % clim([max(climlim)-10 max(climlim)])
        % pbaspect([1 1 1])
       % 
        scatter(AllParams(indsp,6),AllParams(indsp,5),10,log10(valsp),'filled')
        hold on
        try
            scatter(ParamsTest(jChk,6),ParamsTest(jChk,5),100,'ks','filled')
        catch
        end
        scatter(AllParams(inds(end),6),AllParams(inds(end),5),100,'m^','filled')
        IndChosen(jChk)=inds(end);
        xlim([0 0.5])
        ylim([0 0.5])
        if (jChk==nData)
        xlabel('$D_f$')
        else
            xticklabels('')
        end
        ylabel('$k_\textrm{diss}$')
        colormap jet
        climlim=clim;
        clim([max(climlim)-10 max(climlim)])
        pbaspect([1 1 1])
       % 
       %  nexttile
       %  scatter(sqrt(AllParams(indsp,10)),AllParams(indsp,5),10,log10(valsp),'filled')
       %  hold on
       %  try
       %  scatter(sqrt(ParamsTest(jChk,10)),ParamsTest(jChk,5),100,'ks','filled')
       %  catch
       %  end
       %  scatter(sqrt(AllParams(inds(end),10)),AllParams(inds(end),5),100,'m^','filled')
       % %colorbar
       %  climlim=clim;
       %  clim([max(climlim)-10 max(climlim)])
       %  if (jChk==nData)
       %  xlabel('$\ell$')
       %  else
       %      xticklabels('')
       %  end
       %  yticklabels('')
       %  set(gca, 'xScale','Log');
       %  pbaspect([1 1 1])
    end
end