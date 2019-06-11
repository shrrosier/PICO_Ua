function [T0,S0]=PICO_GetOceanInputFromBox0(UserVar,CtrlVar,MUA,GF,ShelfID,PICO_opts)
% Calculate the ocean input for each ice shelf from the ocean input 
% values per basin, weighted with the area of the ice shelf within that
% basin

T0 = zeros(max(ShelfID),1);
S0 = zeros(max(ShelfID),1);

load(PICO_opts.BasinsFile); % Fbasins

if ~strcmp(Fbasins.Method,'nearest')
    Fbasins.Method = 'nearest';
    Fbasins.ExtrapolationMethod = 'nearest';
    warning('Changing basin interpolant to nearest neighbour to ensure unique values');
end

basins = Fbasins(MUA.coordinates(:,1), MUA.coordinates(:,2));

unique_basins = unique(basins);
Nbasins = numel(unique_basins);

MUA2 = MUA;
MUA2.nip = 1;
TriArea = FEintegrate2D([],MUA2,ones(MUA.Nnodes,1));

for ii = 1:max(ShelfID)
    Ak = zeros(Nbasins,1);
    % calc the area of the ice shelf lying in each basin...
    for basin_i = 1:Nbasins
        ind = ShelfID==ii & basins==unique_basins(basin_i);
        Ak(basin_i)=sum(TriArea(ind));
    end
    T0(ii) = dot(PICO_opts.Tbasins, 1/sum(Ak)*Ak);
    S0(ii) = dot(PICO_opts.Sbasins, 1/sum(Ak)*Ak);
end

end