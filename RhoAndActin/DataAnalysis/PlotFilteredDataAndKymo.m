% Code to produce data plots for Fig. 1 in the paper
% It will produce a 2 x 2 plot with the raw data from a movie, followed by
% the filtered data, then the square region and kymograph
% Don't just blindly run this file. Check the threshold for declaring what
% is background
pxlSize=0.1;
figure;
tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact')
NameAll='/Users/ondrejmaxian/Documents/CHICAGO/SRD/Baixue_data/nmy-cyk';
theFiles = dir(fullfile(strcat(NameAll,'/')));
theFolders = theFiles(3:end);
for iF=4%2:length(theFolders)
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
    nexttile
    imagesc(RhoTimeSeries(:,:,1))
    PlotAspect
    colormap sky
    Filtx=smoothdata(RhoTimeSeries,1,'sgolay',50);
    Filtxy=smoothdata(Filtx,2,'sgolay',50);
    FiltRho=smoothdata(Filtxy,3,'sgolay',50);
    nexttile
    imagesc(FiltRho(:,:,1))
    % FiltEm=FiltRho;
    % PlotAspect
    hold on
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
    % Cut out square 200 x 200 region in the center
    xBd1 = find(AllZCols==0,1,'first');
    xBd2 = find(AllZCols==0,1,'last');
    Midx=floor(1/2*(xBd1+xBd2));
    xSt = Midx-WindowSize/2;
    xEnd = Midx+WindowSize/2-1;
    AllZCols(1:xSt-1)=1;
    AllZCols(xEnd+1:end)=1;
    yBd1 = find(AllZRows==0,1,'first');
    yBd2 = find(AllZRows==0,1,'last');
    Midy=floor(1/2*(yBd1+yBd2));
    ySt = Midy-WindowSize/2;
    yEnd = Midy+WindowSize/2-1;
    AllZRows(1:ySt-1)=1;
    AllZRows(yEnd+1:end)=1;
    MeanBoth(AllZRows,:)=0;
    MeanBoth(:,AllZCols)=0;
    Rho=RhoTimeSeries(AllZRows==0,AllZCols==0,:);
    Actin=ActinTimeSeries(AllZRows==0,AllZCols==0,:);
    Info.SaveRows = find(AllZRows==0);
    Info.SaveCols = find(AllZCols==0);
    Info.Name = Name;
    Info.TimeInt= TimeInt;
    Info.IsSD=SD;
    plot(Info.SaveCols(1)*ones(1,200),Info.SaveRows,':k')
    plot(Info.SaveCols(end)*ones(1,200),Info.SaveRows,':k')
    plot(Info.SaveCols,Info.SaveRows(1)*ones(1,200),':k')
    plot(Info.SaveCols,Info.SaveRows(end)*ones(1,200),':k')
    GlobalMeanRho=mean(Rho(:));
    [ny,nx,nFr]=size(Rho);
    for iT=1:nFr
        Rho(:,:,iT)=Rho(:,:,iT)-mean(mean(Rho(:,:,iT)))+GlobalMeanRho;
    end
    Nx=200;
    x=(0:Nx-1)*pxlSize;
    y=(0:Nx-1)*pxlSize;
    % Filter data
    Filtx=smoothdata(Rho,1,'sgolay',50);
    Filtxy=smoothdata(Filtx,2,'sgolay',50);
    FiltRho=smoothdata(Filtxy,3,'sgolay',50);
    nexttile
    imagesc(x,y,FiltRho(:,:,1))
    % Kymograph
    SliceNum=160;
    hold on
    plot([0 nx-1]*pxlSize,pxlSize*SliceNum*[1 1],':k')
    PlotAspect
    nexttile
    RhoKymo=reshape(FiltRho(SliceNum,:,:),200,nFr)';
    imagesc(0.1*(0:199),0.6*(0:nFr),RhoKymo)
    PlotAspect
    pbaspect([1 1 1])
end
