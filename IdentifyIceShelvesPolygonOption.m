function [ShelfID,PBOX,Ak,floating] = IdentifyIceShelvesPolygonOption(CtrlVar,MUA,GF,minArea,minNumS,nmax,MeshBoundaryCoordinates,FloatingCriteria)

%load('Schmitko.mat');
%load('MeshBoundaryCoordinates.mat');
% load ResultsFiles/AntarcticSetup_BFRN_transient_run3/MeshBoundaryCoordinates.mat;

[ShelfID,ShelfGLx,ShelfGLy,ShelfFrontx,ShelfFronty] = GetShelfID(CtrlVar, MUA, GF, MeshBoundaryCoordinates,minArea,minNumS);


Ak = zeros(max(ShelfID),nmax);
rho_m = Ak;

x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);

blnkBox = zeros(MUA.Nnodes,1);


for ii = 1:max(ShelfID)
    
    xShelf = x(ShelfID==ii); yShelf = y(ShelfID==ii);
    
    dGL{ii}=DistanceCloudAtoCloudB(xShelf,yShelf,ShelfGLx{ii},ShelfGLy{ii});
    
    dIF{ii}=DistanceCloudAtoCloudB(xShelf,yShelf,ShelfFrontx{ii},ShelfFronty{ii});
    
%     [daMin,ia,dbMin,ib,D2]=DistanceCloudAtoCloudB(xps,yps,ShelfFrontx{ii},ShelfFronty{ii});
%     T0(ii) = mean(TT(ib));
%     S0(ii) = mean(SS(ib));
    


end

%%

for ii = 1:max(ShelfID)

    
    d_max = max(vertcat(dGL{:}));
    dGL_D = max(dGL{ii});
    nD = 1 + round(sqrt(dGL_D/d_max)*(nmax-1));
    
    rbox = dGL{ii}./(dGL{ii}+dIF{ii});
    
    PBox2 = zeros(numel(dGL{ii}),1);
    
    for k = 1:nD
        
        p1 = 1-sqrt((nD-k+1)/nD);
        p2 = 1-sqrt((nD-k)/nD);
        
        PBox2(p1 <= rbox & rbox <= p2) = k;
    end
    
    
    blnkBox(ShelfID==ii) = PBox2;
    
end

%% THIS SECTION DETERMINES THE BOX AREAS

PBoxEle=round(SNodes2EleMean(MUA.connectivity,blnkBox));
ShelfIDEle = round(SNodes2EleMean(MUA.connectivity,ShelfID));
[Areas,~,~,~]=TriAreaFE(MUA.coordinates,MUA.connectivity);

% Each row of Ak is a unique shelf and each column is a box number within
% that shelf, each element of Ak is the total area of a box in a shelf
for ii = 1:max(ShelfID)
    for k = 1:nmax
       ind = ShelfIDEle==ii & PBoxEle==k;
       Ak(ii,k) = sum(Areas(ind));
    end
end


PBOX = blnkBox;
PBOX(PBOX==0) = nan;

switch FloatingCriteria
    case 'GLthreshold'
    floating = GF.node < CtrlVar.GLthreshold;
    case 'StrictDownstream'
    GF=IceSheetIceShelves(CtrlVar,MUA,GF,[],[],[]);
    floating = GF.NodesDownstreamOfGroundingLines;
    otherwise
    error('Invalid value for PICO_opts.FloatingCriteria');
end

end