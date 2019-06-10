function [ShelfID,ShelfGLx2,ShelfGLy2,ShelfFrontx2,ShelfFronty2] = GetShelfID(CtrlVar, MUA, GF,PICO_opts)
% Identify each ice shelf with an individual ID number
%
% Example:
%
%  [ShelfID,midx,midy,ShelfGLx,ShelfGLy,ShelfFrontx,ShelfFronty] = IdentifyIceShelves_sr(CtrlVar, MUA, GF, MeshBoundaryCoordinates)
%
% - ShelfID: is a vector of length MUA.Nnodes and each ice shelf has a
% unique ID number
%
% current bottlenecks are almost entirely due to potentially larger number
% of calls to inpoly depending on number of seperate GLs and the initial
% call to ArrangeGroundingLinePos which is slow...
%
% possible issues to check... does this properly deal with  (ie. ignore)
% lakes close to the mesh boundary? (within the dist_tol)


% declare BCn as persistent so any reordering only needs to be done once
if PICO_opts.persistentBC
    persistent BCn
end

%==========================================================================
% TUNING PARAMETERS:
GL_cutoff = 0; % if a GL polyline is shorter than this length, ignore it NB: runtime sensitive to this number, 0 will give best behaviour but larger number will make it considerably faster
% minArea = 1e9; % if a shelf area is less than this, ignore it
% minNumShelf = 20; % if a shelf has less floating nodes than this number, ignore it (NOTE: if this is set to zero then grounded islands touching the boundary will be included in the result!!)
%==========================================================================

x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);

[GLgeo,~,~]=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
[xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo);

% create index of where this vector starts a new GL
if any(isnan(xGL))
    GL_ind = [0; find(isnan(xGL))];
else
    GL_ind = [0; numel(xGL)];
end

if PICO_opts.persistentBC && ~isempty(BCn)
    if PICO_opts.InfoLevel > 10
        disp('Using previous boundary coordinate file');
    end
else
    BCn = PICO_opts.MeshBoundaryCoordinates;
end

dxs = x - BCn(1,1); dys = y - BCn(1,2);
ds = hypot(dxs,dys); [~,nod_start] = min(ds); % find index of node at start of MeshBoundaryCoords to check if it's floating
I =  GF.node < CtrlVar.GLthreshold;

if GF.node(nod_start) < 0.5 % now we have a problem - mesh boundary coordinates intersect a shelf
    % the aim here is to find the longest stretch of grounded nodes and put
    % the mesh boundary coordinate intersection here where it will cause no
    % problems... maybe excessive but not too slow
    warning('Mesh boundary coordinates seem to intersect an ice shelf, rearranging start and end points...');
    ib = zeros(length(BCn),1);
    for ii = 1:length(BCn)
        dist =hypot(x-BCn(ii,1),y-BCn(ii,2));
        ib(ii) = find(dist==0);
    end
    GFboundary = GF.node(ib);
    GFchange = diff(GFboundary);
    GFbound_length = 0;
    GFbound_length_max = 0;
    GFbound_pos = nan;
    start_counting = false;
    for ii = 1:numel(GFchange)
        if GFchange(ii) > 0.5 % next node is grounded
            start_counting = true;
            start_ind = ii+1;
        elseif GFchange(ii) < -0.5 % grounded to floating
            start_counting = false;
            stop_ind = ii;
            if GFbound_length > GFbound_length_max
                GFbound_length_max = GFbound_length;
                GFbound_pos = ceil((start_ind+stop_ind)/2);
            end
            GFbound_length = 0;
        end
        if start_counting == true
            GFbound_length = GFbound_length + 1;
        end
    end
    BCn = [BCn(GFbound_pos+1:end,:);BCn(1:GFbound_pos,:)];
end


%initiate blank shelfID vector
ShelfID = zeros(MUA.Nnodes,1)*nan;
shelfnum = 1;

% calculate total area of the domain
% domain_area = polyarea(BCn(:,1),BCn(:,2));

% calculate length (number of nodes) of each GL polyline
GL_length = diff(GL_ind);


num_GL_ignore = 0;
num_lake_ignore = 0;
num_ShelfArea_ignore = 0;
num_ShelfSize_ignore = 0;

% this vector will end up with ones in the indices of any nodes that need
% to be deleted later
to_purge = zeros(MUA.Nnodes,1);

for ii =1:numel(GL_ind)-1 % loop through every GL polyline
    
    % if the number of vertices in the GL polyline is below threshold,
    % ignore and move on (not break just to count number of occurences..)
    if GL_length(ii) < GL_cutoff
        num_GL_ignore = num_GL_ignore + 1;
        continue
    end
    
    % this checks if the GL polyline defines a lake (faster than LakeOrOcean.m)
    % if it is then move on to the next grounding line
    IN = 0;
    if xGL(GL_ind(ii)+1)-xGL(GL_ind(ii+1)-1) == 0
        if yGL(GL_ind(ii)+1)-yGL(GL_ind(ii+1)-1) == 0 % start and end point are the same... ie. this is an enclosed GL polygon... so is it a lake or island?
            % if the GL polygon encloses floating nodes, it must be a lake
            IN = inpoly([x(I) y(I)],[xGL(GL_ind(ii)+1:GL_ind(ii+1)-1) yGL(GL_ind(ii)+1:GL_ind(ii+1)-1)],[],1);
        end
    end
    if sum(IN) > 0 % if lake
        num_lake_ignore = num_lake_ignore + 1;
        continue
    end
    
    % find the nearest node to start and end of mesh boundary coords
    dxs = xGL(GL_ind(ii)+1)-BCn(:,1); dys = yGL(GL_ind(ii)+1)-BCn(:,2);
    ds = hypot(dxs,dys); [~,nod_start] = min(ds);
    dxe = xGL(GL_ind(ii+1)-1)-BCn(:,1); dye = yGL(GL_ind(ii+1)-1)-BCn(:,2);
    de = hypot(dxe,dye); [~,nod_end] = min(de);
    
    % generate xvec, yvec which define the boundary of the a polygon formed
    % from the GL polyline and mesh boundary polyline
    
    if nod_start > nod_end % check ordering of indices
        temp = nod_start;
        nod_start = nod_end;
        nod_end = temp;
        % create polygon from GL and mesh boundary coordinates
        xvec = [xGL(GL_ind(ii)+1:GL_ind(ii+1)-1); BCn(nod_start:nod_end,1)];
        yvec = [yGL(GL_ind(ii)+1:GL_ind(ii+1)-1); BCn(nod_start:nod_end,2)];
    else
        xvec = [xGL(GL_ind(ii)+1:GL_ind(ii+1)-1); flipud(BCn(nod_start:nod_end,1))];
        yvec = [yGL(GL_ind(ii)+1:GL_ind(ii+1)-1); flipud(BCn(nod_start:nod_end,2))];
    end
    
    % find all floating nodes within the polygon defined by xvec and yvec
    in2 = inpoly([x(GF.node<0.5) y(GF.node<0.5)],[xvec,yvec]);
    
    % convert back to full length indices for simplicity
    in = zeros(MUA.Nnodes,1);
    in(GF.node<0.5) = in2;
    in = logical(in);
    
    % if number of nodes within this polygon is below threshold, ignore it
    % and move on to next GL polyline... NB: add these to the purge vector
    % so that these small ice shelves are removed from larger surrounding
    % shelves (happens in nested shelves for complex geometries)
    if sum(in2)<PICO_opts.minNumShelf
        num_ShelfSize_ignore = num_ShelfSize_ignore + 1;
        to_purge = to_purge+in;
        continue
    end
    
    % calculate area of the polygon
    shelf_area = polyarea(xvec,yvec);
    
    % if the area is below threshold, ignore it and move on as above
    if shelf_area < PICO_opts.minArea
        num_ShelfArea_ignore = num_ShelfArea_ignore + 1;
        to_purge = to_purge+in;
        continue
    end
    
    % keep track of the coordinates of the GL and Ice front (used for
    % generating melt boxes later)
    ShelfGLx{shelfnum} = xGL(GL_ind(ii)+1:GL_ind(ii+1)-1);
    ShelfGLy{shelfnum} = yGL(GL_ind(ii)+1:GL_ind(ii+1)-1);
    ShelfFrontx{shelfnum} = BCn(nod_start:nod_end,1);
    ShelfFronty{shelfnum} = BCn(nod_start:nod_end,2);
    
    % save main results of this polygon which will be reordered later...
    insav{shelfnum} = in;
    midx(shelfnum) = mean(xvec);
    midy(shelfnum) = mean(yvec);
    area_sav(shelfnum) = shelf_area;
    xvecsav{shelfnum} = xvec;
    yvecsav{shelfnum} = yvec;
    
    shelfnum = shelfnum+1; % we got this far so it was a valid ice shelf, increment shelf number
    
end

if shelfnum==1
    error('No valid ice shelves detected - check shelf size and shelf area cutoffs are sensible for your domain');
end

% this section sorts all the identified ice shelves by area (not
% necessarily already the case - previously sorted by number of GL
% vertices). We need to do this to help with the check for nested ice
% shelves.

[~,ia] = sort(area_sav,'descend');
insav2 = insav(ia);
ShelfFrontx2 = ShelfFrontx(ia);
ShelfFronty2 = ShelfFronty(ia);
ShelfGLx2 = ShelfGLx(ia);
ShelfGLy2 = ShelfGLy(ia);
midx2 = midx(ia);
midy2 = midy(ia);

% this section assigns ice shelf IDs to nodes and does more checks for
% nested shelves

for ii = 1:max(shelfnum)-1
    
    ShelfID(insav2{ii}) = ii; % assign shelf IDs
    
    for jj = 1:ii-1 % go back through all previously defined shelves, NB: these are guaranteed to be larger than the one we just created
        if sum(insav2{ii} & insav2{jj}) > 0 % if a newly created shelf shares nodes with a previous shelf, this new shelf must be contained within the larger shelf
            big_shelf = [ShelfFrontx2{jj} ShelfFronty2{jj}]; % ice front coordinates of the surrounding shelf
            small_shelf = [ShelfFrontx2{ii} ShelfFronty2{ii}]; % ice front coordinates of the nested shelf
            [~,not_in_small] = setdiff(big_shelf,small_shelf,'rows'); % coords in big shelf that are not in small shelf
            ShelfFrontx2{jj} = ShelfFrontx2{jj}(not_in_small); % remove any ice front x coords from the larger shelf which were also in the smaller shelf
            ShelfFronty2{jj} = ShelfFronty2{jj}(not_in_small); % same for y coords
            ShelfID(insav2{jj}(not_in_small)) = nan; % also remove these nodes from the previously defined larger ice shelf
        end
    end
    
end

% now we can delete all those nodes belonging to small shelves that didn't
% make the cutoff... I think this cannot be done before this point...
ShelfID(to_purge>0) = nan;

% finally, with all other checks done, we count the number of nodes with
% each unique shelf ID and if the total is below the threshold, delete this
% shelf and renumber shelves with a higher shelf ID
badShelf = [];
cur = 1;
for ii = 1:max(ShelfID)
    numShelf = sum(ShelfID==cur);
    if numShelf<PICO_opts.minNumShelf
        badShelf = [badShelf ii];
        ShelfID(ShelfID==cur) = nan;
        ShelfID(ShelfID>cur) = ShelfID(ShelfID>cur)-1;
    else
        cur = cur+1;
    end
end
midx2(badShelf) = [];
midy2(badShelf) = [];
ShelfGLx2(badShelf) = [];
ShelfGLy2(badShelf) = [];
ShelfFrontx2(badShelf) = [];
ShelfFronty2(badShelf) = [];

if PICO_opts.InfoLevel>10
    
    disp(['Ignored ' num2str(num_GL_ignore) ' ice shelves due to GL length cutoff criteria of ' num2str(GL_cutoff)]);
    disp(['Ignored ' num2str(num_ShelfArea_ignore) ' ice shelves due to shelf area cutoff criteria of ' num2str(PICO_opts.minArea)]);
    disp(['Ignored ' num2str(num_ShelfSize_ignore) ' ice shelves due to shelf size cutoff criteria of ' num2str(PICO_opts.minNumShelf)]);
    disp(['Ignored ' num2str(num_lake_ignore) ' lakes (floating nodes fully enclosed by a GL)']);
    
end
