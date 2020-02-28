function Continent_nodes = get_continental_nodes(CtrlVar,MUA,GF,Areas,PICO_opts)

if ~isfield(GF,'NodesCrossingGroundingLines')
    GF = IceSheetIceShelves(CtrlVar,MUA,GF);
end

EleSubset = GF.ElementsUpstreamOfGroundingLines;

grounded = GF.node>0.5;

TRI=MUA.connectivity(EleSubset,:) ;

% M = connectivity2adjacency(TRI)>0 ;
% G=graph(M) ;

G=graph(TRI,TRI(:,[2 3 1]));
% bins=conncomp(G) ;

Continent_nodes = zeros(MUA.Nnodes,1);
% [Areas,~,~,~]=TriAreaFE(MUA.coordinates,MUA.connectivity);

[bins,binsize]=conncomp(G) ;
badbins = binsize<100;
gd = find(~badbins);
Lia = ismember(bins,gd)';
grounded(~Lia) = 0; % get rid of grounded nodes that are not connected to sufficient neighbours


while sum(grounded)>0
    NodeSeed = find(grounded,1,'first');
    ID=bins(NodeSeed) ;

    nodes=find(bins==ID);
    if numel(nodes) > 100
        blnk = zeros(MUA.Nnodes,1); blnk(nodes) = 1;
        Ele = logical(ceil(Nodes2EleMean(MUA.connectivity,blnk)));
        cont_area = sum(Areas(Ele));
        if cont_area > PICO_opts.ContinentArea
            Continent_nodes(nodes) = 1;
        end
    end
    grounded(nodes) = 0;
end

