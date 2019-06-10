function [ShelfID,PBOX,Ak,floating] = IdentifyIceShelvesPolygonOption(CtrlVar,MUA,GF,PICO_opts)

nmax = PICO_opts.nmax;

if PICO_opts.InfoLevel>10
    disp('===========================================================');
    fprintf('Calling GetShelfID...\n');
end

[ShelfID,ShelfGLx,ShelfGLy,ShelfFrontx,ShelfFronty] = GetShelfID(CtrlVar, MUA, GF, PICO_opts);
if PICO_opts.InfoLevel>10
    fprintf('Polygon option ran successfully and defined %2i ice shelves\n',max(ShelfID));
end
Ak = zeros(max(ShelfID),nmax);

x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);

blnkBox = zeros(MUA.Nnodes,1);

if PICO_opts.InfoLevel>10
    fprintf('For each Ice Shelf node, calculating distance to GL and ice front...');
end

for ii = 1:max(ShelfID)
    
    xShelf = x(ShelfID==ii); yShelf = y(ShelfID==ii);
    
    dGL{ii}=DistanceCloudAtoCloudB(xShelf,yShelf,ShelfGLx{ii},ShelfGLy{ii});
    
    dIF{ii}=DistanceCloudAtoCloudB(xShelf,yShelf,ShelfFrontx{ii},ShelfFronty{ii});
    
%     [daMin,ia,dbMin,ib,D2]=DistanceCloudAtoCloudB(xps,yps,ShelfFrontx{ii},ShelfFronty{ii});
%     T0(ii) = mean(TT(ib));
%     S0(ii) = mean(SS(ib));
    
end

if PICO_opts.InfoLevel>10
    fprintf('Success!\n');
end

%%

if PICO_opts.InfoLevel>10
    fprintf('Separating each ice shelf into a maximum of %1i boxes ',PICO_opts.nmax);
    fprintf('(You can change this value in PICO_opts.nmax)...');
end

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

if PICO_opts.InfoLevel>10
    fprintf('Success!\n');
end

%% THIS SECTION DETERMINES THE BOX AREAS

PBoxEle=round(SNodes2EleMean(MUA.connectivity,blnkBox));
ShelfIDEle = round(SNodes2EleMean(MUA.connectivity,ShelfID));
[Areas,~,~,~]=TriAreaFE(MUA.coordinates,MUA.connectivity);

if PICO_opts.InfoLevel>10
    fprintf('Calculating the area of each box...');
end

% Each row of Ak is a unique shelf and each column is a box number within
% that shelf, each element of Ak is the total area of a box in a shelf
for ii = 1:max(ShelfID)
    for k = 1:nmax
       ind = ShelfIDEle==ii & PBoxEle==k;
       Ak(ii,k) = sum(Areas(ind));
    end
end

if PICO_opts.InfoLevel>10
    fprintf('Success!\n');
end


PBOX = blnkBox;
PBOX(PBOX==0) = nan;

switch PICO_opts.FloatingCriteria
    case 'GLthreshold'
    floating = GF.node < CtrlVar.GLthreshold;
    case 'StrictDownstream'
    GF=IceSheetIceShelves(CtrlVar,MUA,GF,[],[],[]);
    floating = GF.NodesDownstreamOfGroundingLines;
    otherwise
    error('Invalid value for PICO_opts.FloatingCriteria');
end

if PICO_opts.InfoLevel>10
    disp('===========================================================');
    disp('======== Polygon algorithm completes successfully =========');
    disp('===========================================================');
end

end