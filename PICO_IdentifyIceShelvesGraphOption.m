function [ShelfNum,BoxID,Ak,floating] = PICO_IdentifyIceShelvesGraphOption(UserVar,CtrlVar,MUA,GF,PICO_opts)

x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);

GF = IceSheetIceShelves(CtrlVar,MUA,GF);
EleSubset = GF.ElementsDownstreamOfGroundingLines;

[~,OceanNodes] = LakeOrOcean3(CtrlVar,MUA,GF);
floating = OceanNodes;

FloatNodesTemp = floating;

TRI=MUA.connectivity(EleSubset,:) ;

G=graph(TRI,TRI(:,[2 3 1]),TRI*0+1);
% G = simplify(G);

[bins,binsize]=conncomp(G) ;
badbins = binsize<PICO_opts.minNumShelf;
gd = find(~badbins);
Lia = ismember(bins,gd)';
FloatNodesTemp(~Lia) = 0; % get rid of floating nodes that are not connected to sufficient neighbours

ii=1;

ShelfNum = zeros(MUA.Nnodes,1);
while sum(FloatNodesTemp)>0
    NodeSeed = find(FloatNodesTemp,1,'first');
    ID=bins(NodeSeed) ;
    nodes=find(bins==ID);
    %shelf{ii} = nodes;
    ShelfNum(nodes) = ii;
    ii = ii+1;
    FloatNodesTemp(nodes) = 0;
end


% this section could be used to calculate distances but it seems to be very
% slow
%
% FloatingBoundaryNodes=MUA.Boundary.Nodes(GF.node(MUA.Boundary.Nodes)<0.5);
% node_vec = 1:MUA.Nnodes;
% dIF = node_vec*0;
% dGL = dIF;
% 
% 
% for ii = 1:max(ShelfNum)
%     nodes = node_vec(ShelfNum==ii);
%     calv_nodes = intersect(nodes,FloatingBoundaryNodes);
%     main_shelf_nodes = setdiff(nodes,calv_nodes);
%     d = distances(G,main_shelf_nodes,calv_nodes);
%     mind = min(d,[],2);
%     dIF(main_shelf_nodes) = mind;
%     dIF(calv_nodes) = 0;
% end
% 

%%
[Areas,~,~,~]=TriAreaFE(MUA.coordinates,MUA.connectivity);
Continent_nodes1 = get_continental_nodes(CtrlVar,MUA,GF,Areas,PICO_opts);
NodeNums = 1:MUA.Nnodes';
NodeNums(~Continent_nodes1 | ~GF.NodesCrossingGroundingLines)= [];
NodeInEle = ismember(MUA.TR,NodeNums);
EleInd = max(NodeInEle,[],2);
Continent_nodes2 = unique(MUA.TR(EleInd,:));

FloatingBoundaryNodes=MUA.Boundary.Nodes(GF.node(MUA.Boundary.Nodes)<0.5);
dIF = zeros(numel(ShelfNum),1);
dGL = dIF;
NodeNums = 1:MUA.Nnodes';
NodeNums2 = zeros(MUA.Nnodes,1);
NodeNums2(Continent_nodes2) = 1; 


for ii = 1:max(ShelfNum) % calculate  dIF for each ice shelf
    xBox = x(ShelfNum==ii);
    yBox = y(ShelfNum==ii);
    shelf_boundary = ismember(NodeNums(ShelfNum==ii),FloatingBoundaryNodes);
    %[~, D] = knnsearch([xBox(shelf_boundary) yBox(shelf_boundary)],[xBox yBox]); % distance of every shelf node to calving front
    [D,I] = pdist2([xBox(shelf_boundary) yBox(shelf_boundary)],[xBox yBox],'euclidean','Smallest',1); % faster option than knnsearch
    dIF(ShelfNum==ii) = D;
    GLnodes = ShelfNum==ii & NodeNums2;
    if sum(GLnodes)==0
        dGL(ShelfNum==ii) = 1;
        warning('No grounded nodes detected for Ice Shelf %03i ',ii);
        continue
    end
    [~, D2] = knnsearch([x(GLnodes) y(GLnodes)],[xBox yBox]); 
    dGL(ShelfNum==ii) = D2; 
end


ShelfNum(~floating) = nan;

%%

dmax = max(dGL);
BoxID = zeros(size(dGL));
badBox = zeros(max(ShelfNum),1);

for ii = 1:max(ShelfNum) %need a second loop because only now do we know dmax
    ind = ShelfNum==ii;
    dglmax = max(dGL(ind));
    nD = 1 + round(sqrt(dglmax./dmax)*(PICO_opts.nmax-1));
    
    rbox = dGL(ind)./(dGL(ind)+dIF(ind));
    blnkBox = rbox*0;
    
    for k = 1:nD
        p1 = 1-sqrt((nD-k+1)/nD);
        p2 = 1-sqrt((nD-k)/nD);
        blnkBox(p1 <= rbox & rbox <= p2) = k;
        if k == 1 && sum(blnkBox==1) ==0 % this checks if any ice shelves are missing box #1 which might occasionally happen for various reasons
            warning('Ice Shelf %03i is missing a box ID # 1, this might be because of problems with PICO options, in particular check the continent cutoff value',ii);
            badBox(ii) = 1;
        end
    end
    
    BoxID(ind) = blnkBox;
    
end

if any(badBox==1)
    badInds = find(badBox==1);
    for ii = 1:numel(badInds)
        ind = ShelfNum==badInds(ii) & BoxID == 1;
         while sum(ind)==0
            BoxID(ShelfNum==badInds(ii)) = BoxID(ShelfNum==badInds(ii))-1;
            ind = ShelfNum==badInds(ii) & BoxID == 1;
            if sum(ind)>0
                badBox(badInds(ii)) = 0;
            end
         end
    end
end

Ak = PICO_calc_areas(MUA, BoxID, ShelfNum, PICO_opts.nmax);
