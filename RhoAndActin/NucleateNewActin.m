% Nucleate new actin 
function AddedNucs=NucleateNewActin(NucRate,x,y)
    % Nucleate a filament 
    dx=x(2)-x(1);
    Nx=length(x);
    dy=y(2)-y(1);
    Ny=length(y);
    NumNuc = NucRate*dx*dy;
    nDef = floor(NumNuc);
    nAdd = rand(Ny,Nx) < (NumNuc-nDef);
    nTot = nDef+nAdd; % number of nucleates for each cell
    NewCells = find(nTot>0);
    AddedNucs=[];
%     nMax=max(nTot(:));
%     AddedNucs2=zeros(length(NewCells)*nMax,2)+nan;
    for iC=1:length(NewCells)
        yind=mod(NewCells(iC),Ny);
        if (yind==0)
            yind=Ny;
        end
        xind = ceil(NewCells(iC)/Ny-1e-10);
        nForCell = nTot(yind,xind);
        if (nForCell ==0)
            keyboard
        end
        % Nucleate in the box between x and x+dx
        NewPt = rand(1,2).*[dx dy]+[x(xind) y(yind)];
        AddedNucs=[AddedNucs;NewPt];
%         AddedNucs2((iC-1)*nMax+1:(iC-1)*nMax+nTot(NewCells(iC)),:)=NewPt;
%         try
%         AddedNucs2=rmmissing(AddedNucs2);
%         if(max(abs(AddedNucs2(:)-AddedNucs(:)))>0)
%             keyboard
%         end
%         catch
%             keyboard
%         end
    end
end
