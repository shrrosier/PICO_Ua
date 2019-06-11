function [ShelfNum,BoxID,Ak,floating] = PICO_MakeBoxesOneShelf(UserVar,CtrlVar,MUA,GF,PICO_opts)

%% first create the ice shelf ID from a floating mask

nmax = PICO_opts.nmax;
FloatingCriteria = PICO_opts.FloatingCriteria;

switch FloatingCriteria
    case 'GLthreshold'
        floating = GF.node < CtrlVar.GLthreshold;
    case 'StrictDownstream'
        GF=IceSheetIceShelves(CtrlVar,MUA,GF,[],[],[]);
        floating = GF.NodesDownstreamOfGroundingLines;
    otherwise
        error('Invalid value for PICO_opts.FloatingCriteria');
end
ShelfNum = floating;


%% this section calculates distances from GLs and ice fronts for each shelf
FloatingBoundaryNodes=MUA.Boundary.Nodes(GF.node(MUA.Boundary.Nodes)<0.5);
dIF = ShelfNum*0;
dGL = dIF;

x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);

xBox = x(ShelfNum==1);
yBox = y(ShelfNum==1);
[~, D] = knnsearch([x(FloatingBoundaryNodes) y(FloatingBoundaryNodes)],[xBox yBox]); % distance of every shelf node to calving front
dIF(ShelfNum==1) = D;

% Calculate distance from GL... this potentially will slow things a
% lot

% option 1: use boundary to get bounding nodes for each ice shelf and
% then exclude those nodes that are also part of MUA.Boundary.Nodes
% then use knnsearch on this much smaller subset of nodes to get
% distance
% advantage of removing islands within an ice shelf

k = boundary(xBox,yBox,1);
te=setdiff(k,find(D==0)); % instead of this... ensure node is not part of the same element as a node on the boundary
[~, D2] = knnsearch([xBox(te) yBox(te)],[xBox yBox]); % distance of every shelf node to GL
dGL(ShelfNum==1) = D2;

% option 2: just search over all GL nodes
% this is actually faster but need a way to get rid of islands from
% the search
%     [~, D2] = knnsearch([x(GLnodes) y(GLnodes)],[xBox yBox]); % distance of every shelf node to calving front
%     dGL(PBOX==ii) = D2;




%% now calculate the box numbers for each ice shelf

dmax = max(dGL);
BoxID = zeros(size(dGL));

    ind = ShelfNum==1;
    dglmax = max(dGL(ind));
    nD = 1 + round(sqrt(dglmax./dmax)*(nmax-1));
    
    rbox = dGL(ind)./(dGL(ind)+dIF(ind));
    blnkBox = rbox*0;
    
    for k = 1:nD
        p1 = 1-sqrt((nD-k+1)/nD);
        p2 = 1-sqrt((nD-k)/nD);
        blnkBox(p1 <= rbox & rbox <= p2) = k;
    end
    
    BoxID(ind) = blnkBox;


%% finally calculate the area of each box in each ice shelf

Ak = zeros(max(ShelfNum),nmax);

PBoxEle=ceil(SNodes2EleMean(MUA.connectivity,BoxID));
ShelfIDEle = round(SNodes2EleMean(MUA.connectivity,ShelfNum));
[Areas,~,~,~]=TriAreaFE(MUA.coordinates,MUA.connectivity);

% Each row of Ak is a unique shelf and each column is a box number within
% that shelf, each element of Ak is the total area of a box in a shelf
    for k = 1:nmax
        ind = ShelfIDEle==1 & PBoxEle==k;
        Ak(1,k) = sum(Areas(ind));
    end
