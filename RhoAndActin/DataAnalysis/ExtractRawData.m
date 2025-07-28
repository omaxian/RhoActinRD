% Extract a 20 x 20 um (200 pxl x 200 pxl) box in the
% center
NameAll='nmy-cyk';
theFiles = dir(fullfile(strcat(NameAll,'/')));
theFolders = theFiles(3:end);
for iF=2:length(theFolders)
    % Go into the folder
    Name = theFolders(iF).name;
    RhoFile = dir(fullfile(strcat(NameAll,'/',Name),strcat(Name,'*','G.tif')));
    ActinFile = dir(fullfile(strcat(NameAll,'/',Name),strcat(Name,'*','R.tif')));
    SD=0;
    if (contains(RhoFile.name,"SD"))
        SD=1;
        TimeInt=0.6;
    else
        TimeInt=0.4; %TIRF
    end
    WindowSize=200;
    RhoTimeSeries=tiffreadVolume(RhoFile.name);
    ActinTimeSeries=tiffreadVolume(ActinFile.name);
    [nY,nX,nFrR]=size(RhoTimeSeries);
    [~,~,nFrA]=size(ActinTimeSeries);
    minFr=min(nFrR,nFrA);
    RhoTimeSeries=RhoTimeSeries(:,:,1:minFr);
    ActinTimeSeries=ActinTimeSeries(:,:,1:minFr);
    % Pull square region in center of embryo
    MeanRho = mean(RhoTimeSeries,3);
    MeanActin =mean(ActinTimeSeries,3);
    MeanBoth = 1/2*(MeanRho+MeanActin);
    MeanBothOG=MeanBoth;
    % Cut out very dim parts
    Bkgrnd=min(MeanBoth(:));
    MaxInt=max(MeanBoth(:));
    % This is the threshold part (VERY IMPORTANT to get good data).
    % For TIRF set to 0.1, for SD 0.4
    if (SD)
        BgThres = 0.4;
    else
        BgThres = 0.1;
    end
    MeanBoth(MeanBoth < Bkgrnd+BgThres*(MaxInt-Bkgrnd))=0;
    % Cut out rows and cols that are all 0 (center embryo)
    AllZCols = sum(MeanBoth==0)==nY;
    AllZRows = sum(MeanBoth'==0)==nX;
    tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
    nexttile
    imagesc(MeanBoth)
    % Cut out square 200 x 200 region in the center
    xBd1 = find(AllZCols==0,1,'first');
    xBd2 = find(AllZCols==0,1,'last');
    Midx=floor(1/2*(xBd1+xBd2));
    xSt = Midx-WindowSize/2;
    xEnd = Midx+WindowSize/2-1;
    if (xSt < xBd1 || xEnd > xBd2)
        keyboard
    end
    AllZCols(1:xSt-1)=1;
    AllZCols(xEnd+1:end)=1;
    yBd1 = find(AllZRows==0,1,'first');
    yBd2 = find(AllZRows==0,1,'last');
    Midy=floor(1/2*(yBd1+yBd2));
    ySt = Midy-WindowSize/2;
    yEnd = Midy+WindowSize/2-1;
    if (ySt < yBd1 || yEnd > yBd2)
        keyboard
    end
    AllZRows(1:ySt-1)=1;
    AllZRows(yEnd+1:end)=1;
    MeanBoth(AllZRows,:)=0;
    MeanBoth(:,AllZCols)=0;
    nexttile
    imagesc(MeanBoth)
    drawnow
    RhoData=RhoTimeSeries(AllZRows==0,AllZCols==0,:);
    ActinData=ActinTimeSeries(AllZRows==0,AllZCols==0,:);
    Info.SaveRows = find(AllZRows==0);
    Info.SaveCols = find(AllZCols==0);
    Info.Name = Name;
    Info.TimeInt= TimeInt;
    Info.IsSD=SD;
    %save(strcat(NameAll,'Info_',num2str(iF)),'Info');
    %save(strcat(NameAll,'Actin_',num2str(iF)),'ActinData');
    %save(strcat(NameAll,'Rho_',num2str(iF)),'RhoData');
end
