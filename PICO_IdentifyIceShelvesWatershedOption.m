function [ShelfNum,BoxID,Ak,floating] = PICO_IdentifyIceShelvesWatershedOption(UserVar,CtrlVar,MUA,GF,PICO_opts)
%
% Function to generate unique shelf IDs with corresponding areas,
% subdivided into boxes using the method described in Reese (2018).
%
% Usage: [ShelfNum,BoxID,ShelfArea] = IdentifyIceShelvesWatershedOption(CtrlVar,MUA,GF,PICO_opts.PICOres,PICO_opts.minArea,minNumS,nmax)
%
% PICO_opts.PICOres = 10000; % the resolution used in the watershed algorithm, lower numbers will be slower but will capture more detailed ice shelf geometries that might be missed otherwise
% PICO_opts.minArea = 2e9; % minimum ice shelf area in m^2 - lower numbers may slow code quite significantly
% minNumS = 20; % minimum number of floating nodes needed to be considered
% an ice shelf, this is only needed to deal with a few cases where an ice shelf may cover a large area in the structured grid but only actually consist of a few nodes
% nmax= 5; % maximum number of boxes per ice shelf


%% this section calculates a unique ID for each shelf using the watershed function

% PICOres = 10000; % the resolution used in the watershed algorithm, will potentially slow things down a lot if this is increased
% PICO_opts.minArea = 1e9; % minimum shelf area in m^2
% nmax= 5; % maximum number of boxes per ice shelf

[~,~,Zi,in] = tri2grid(MUA,GF.node,PICO_opts.PICOres);

notshelf = Zi>0.99 | ~in;
% ws = watershed(notshelf); %bwlabel seems to do a better and MUCH faster job
inshelf = ~notshelf;
ws = bwlabel(inshelf,8);
% CC = bwconncomp(inshelf,8); ws = labelmatrix(CC);  this seems to be slower

grounded_regions = bwlabel(Zi>0.99,8);
gr_vec = reshape(grounded_regions,[],1);
[num_occur,ind] = hist(gr_vec,unique(gr_vec));

continent_cutoff = PICO_opts.ContinentArea/PICO_opts.PICOres^2;

gdGL = num_occur>continent_cutoff & ind'>0; %ignore background (b=0)
gdGLind = ind(gdGL);

loc1 = ismember(grounded_regions,gdGLind); %matrix where 1 = 'continental' grounded ice and 0 is islands/ocean/ice shelf

dGLmat = bwdist(loc1);
x = MUA.coordinates(:,1); y= MUA.coordinates(:,2);

%x2 and y2 are coordinates of bottom left corner of pixels
x2 = floor((x-min(x))./PICO_opts.PICOres) + 1;
y2 = floor((y-min(y))./PICO_opts.PICOres) + 1;

ShelfID = zeros(MUA.Nnodes,1);
dGL = ShelfID;

for ii = 1:numel(x)
    ShelfID(ii) = max([ws(y2(ii),x2(ii)) ws(y2(ii)+1,x2(ii)) ws(y2(ii),x2(ii)+1) ws(y2(ii)+1,x2(ii)+1)]);
    dGL(ii) = min([dGLmat(y2(ii),x2(ii)) dGLmat(y2(ii)+1,x2(ii)) dGLmat(y2(ii),x2(ii)+1) dGLmat(y2(ii)+1,x2(ii)+1)]) .* PICO_opts.PICOres;
end

% possibly add here - if any triangle has a node beloning to an ice shelf
% then make all nodes of that triangle belong to the same ice shelf (this
% will hopefully deal with some edge issues

switch PICO_opts.FloatingCriteria
    case 'GLthreshold'
    floating = GF.node < CtrlVar.GLthreshold;
    case 'StrictDownstream'
    GF=IceSheetIceShelves(CtrlVar,MUA,GF,[],[],[]);
    floating = GF.NodesDownstreamOfGroundingLines;
    otherwise
    error('Invalid value for PICO_opts.FloatingCriteria');
end
ShelfID(~floating) = -1;

%% this section roughly calculates the area of each ice shelf and removes ice shelves below some minimum area
ws2 = ws; ws2(notshelf)=nan; ws2 = reshape(ws2,[],1);
numshelf = histc(ws2,1:max(max(ws2)));
numshelf2 = histc(ShelfID,1:max(ShelfID));

minNumG = round(PICO_opts.minArea/PICO_opts.PICOres^2);
badShelf1 = find(numshelf<minNumG);
badShelf2 = find(numshelf2<PICO_opts.minNumShelf);
badShelf = union(badShelf1,badShelf2);
idx = ismember(ShelfID,badShelf); % find the indices of shelves considered bad
ShelfID(idx) = -1; % replace the shelf ID in these locations
ShelfID(ShelfID==0) = -1; % watershed algorithm gives 0 along watershed boundaries, not sure of the best way to deal with this right now...
[~,~,unShelfID] = unique(ShelfID); % replaces all shelf IDs with ascending numbers from 1 (note at this point 1 is grounded ice)
unShelfID(unShelfID<2) = nan;
ShelfNum = unShelfID-1;


if max(ShelfID)==-1
    error('No valid ice shelves detected - check shelf size and shelf area cutoffs are sensible for your domain');
end


%% this section calculates distances from ice fronts for each shelf
FloatingBoundaryNodes=MUA.Boundary.Nodes(GF.node(MUA.Boundary.Nodes)<0.5);
dIF = zeros(numel(ShelfNum),1);

for ii = 1:max(ShelfNum) % calculate  dIF for each ice shelf
    xBox = x(ShelfNum==ii);
    yBox = y(ShelfNum==ii);
    [~, D] = knnsearch([x(FloatingBoundaryNodes) y(FloatingBoundaryNodes)],[xBox yBox]); % distance of every shelf node to calving front
    dIF(ShelfNum==ii) = D;
end

%% now calculate the box numbers for each ice shelf

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

% if any shelves don't have a box # 1, this loop renumbers the boxes so
% until there is at least one node in box # 1
%
% I don't think it's possible for an ice shelf to skip a box number ie have
% a box # 1 and # 3 but no box #2
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

%% finally calculate the area of each box in each ice shelf

Ak = PICO_calc_areas(MUA, BoxID, ShelfNum, PICO_opts.nmax);

end