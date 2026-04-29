% Process and filter frog data
MasterFolder = strcat('/Users/ondrejmaxian/Documents/SRD/',...
    'Frog data for Ondrej/Frog stuff for Ondrej/');
Doses = [33 66 166 333 1000];
for jD=1:length(Doses)
    TIFList{jD}=[];
    DoseFolder = strcat(num2str(Doses(jD)),' ng copy/');
    subfold=dir(strcat(MasterFolder,DoseFolder));
    for sF=1:length(subfold)
        if (contains(subfold(sF).name,'6s'))
            % Get all the tifs in that folder
            insidesub = dir(strcat(MasterFolder,DoseFolder,subfold(sF).name));
            for k=1:length(insidesub)
                if (contains(insidesub(k).name,'.tif') &&~contains(insidesub(k).name,'BAD'))%  && contains(insidesub(k).name,'_A-001'))
                    TIFList{jD}=[TIFList{jD};convertCharsToStrings(insidesub(k).name)];
                end
            end
        end
    end
end

nDose = length(Doses);
RepIndexByDose = [7 2 16 1 10];
rng(0);
tiledlayout(3,length(Doses),'Padding', 'none', 'TileSpacing', 'compact');
nModesFilt=100;
FrTime = 6;
pxlSize = 0.2661449; % um
AllNames=strings(length(Doses),1);
for jD=1:nDose
    % Pick out a random movie
    ListForDose = TIFList{jD};
    if ~isempty(ListForDose)
    %figure;
    end
    %index = ceil(rand*length(ListForDose));
    for index=RepIndexByDose(jD)
    MovieName = convertStringsToChars(ListForDose(index));
    AllNames(jD)=convertCharsToStrings(MovieName);
    Comb=tiffreadVolume(MovieName);
    Rho=Comb(:,:,1:2:end);
    Actin=Comb(:,:,2:2:end);
    % % Difference subtraction (removes static signal)
    MF = floor(size(Rho,3)/2);
    Rho=Rho-Rho(:,:,end);
    Actin=Actin-Actin(:,:,end);
    Rho=double(Rho(:,:,1:end-1));
    Actin=double(Actin(:,:,1:end-1));
    [FiltRho,RhoHatMean,ACorsRho,TimeLags] = FilterData(Rho,nModesFilt,FrTime);
    [FiltActin,ActHatMean,ACorsAct,TimeLags] = FilterData(Actin,nModesFilt,FrTime);
    % Find the % of actin that colocalizes with Rho
    Spots=FiltRho>prctile(FiltRho(:),90);
    prcColoc = sum(FiltActin(Spots(:)))/sum(FiltActin(:));
    Nx = size(Rho,2);
    Ny = size(Rho,1);
    x=(0:Nx-1)*pxlSize;
    y=(0:Ny-1)*pxlSize;
    % Compute XCor
    [UvalsF,dtvalsF,DistsByRF] = CrossCorrelations(pxlSize,pxlSize,FrTime,...
        FiltRho,FiltActin,1);
    DistsByRF=DistsByRF/max(abs(DistsByRF(:)));
    ResampledT = -120:2:120;
    ResampledX = 0:0.5:10;
    [X,T]=meshgrid(ResampledX,ResampledT);
    X0Inds=find(X==0);
    T0Inds=find(T==0);
    DataXCor=ResampleXCor(DistsByRF,dtvalsF,UvalsF,ResampledX,ResampledT,11,121);
    % Make a kymograph
    %figure
    %tiledlayout(3,1,'Padding', 'none', 'TileSpacing', 'compact');
    ax0=nexttile(jD);
    imagesc((0:Nx-1)*pxlSize,(0:Nx-1)*pxlSize,FiltRho(:,:,MF))
    clim([min(FiltRho(:)) max(FiltRho(:))])
    colormap(ax0,sky)
    pbaspect([1 1 1])
    %xlabel('$x$ ($\mu$m)')
    ax1=nexttile(jD+nDose);
    C2=[0.87 0.49 0];
    C1=[0.95 0.9 0.9];
    Cmap=C1+(0:100)'/100.*(C2-C1);
    RhoKymo=reshape(FiltRho(2*Ny/4,:,:),Nx,size(FiltRho,3))';
    ActKymo=reshape(FiltActin(2*Ny/4,:,:),Nx,size(FiltActin,3))';
    imagesc((0:Nx-1)*pxlSize,(0:size(RhoKymo,1)-1)*FrTime,RhoKymo)
    clim([min(FiltRho(:)) max(FiltRho(:))])
    %xlabel('$x$ ($\mu$m)')
    %ylabel('$t$ (s)')
    pbaspect([1 1 1])
    colormap(ax1,sky)
    ax2=nexttile(jD+2*nDose);
    imagesc(ResampledX,ResampledT,DataXCor)
    clim([-1 1])
    %xlabel('$\Delta r$ ($\mu$m)')
    %ylabel('$\Delta t$ (s)')
    pbaspect([1 1 1])
    colormap(ax2,turbo)
    %Full movie 
    % figure(2*jD)
    % for g=1:size(Rho,3)
    %     tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
    %     ax3=nexttile;
    %     imagesc(x,y,FiltRho(:,:,g))
    %     pbaspect([1 1 1])
    %     colormap(ax3, sky)
    %     clim([min(FiltRho(:)) max(FiltRho(:))])
    %     ax4=nexttile;
    %     imagesc(x,y,FiltActin(:,:,g))
    %     pbaspect([1 1 1])
    %     colormap(ax4,Cmap)
    %     clim([min(FiltActin(:)) max(FiltActin(:))])
    %     drawnow
    % end
    end
    XCorsByDose(:,:,jD)=DataXCor;
end
