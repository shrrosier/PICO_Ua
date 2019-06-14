function [T0,S0]=PICO_GetOceanInputFromBox0(UserVar,CtrlVar,MUA,GF,ShelfID,PICO_opts)
% Calculate the ocean input for each ice shelf from the ocean input
% values per basin, weighted with the area of the ice shelf within that
% basin

T0 = zeros(max(ShelfID),1);
S0 = zeros(max(ShelfID),1);

if ~isfield(PICO_opts,'BasinsInterpolant')
    basins = ones(MUA.Nnodes,1);
else
    if ~strcmp(PICO_opts.BasinsInterpolant.Method,'nearest')
        PICO_opts.BasinsInterpolant.Method = 'nearest';
        PICO_opts.BasinsInterpolant.ExtrapolationMethod = 'nearest';
        warning('Changing basin interpolant to nearest neighbour to ensure unique values');
    end
    basins = PICO_opts.BasinsInterpolant(MUA.coordinates(:,1), MUA.coordinates(:,2));
end

unique_basins = unique(basins);
Nbasins = numel(unique_basins);

if ~isfield(PICO_opts,'Tbasins')
    defaultT = -1.8; % deg C
    PICO_opts.Tbasins = zeros(Nbasins,1)+defaultT;
    if PICO_opts.InfoLevel>0
        warning(strcat('Ocean input vector is missing, setting temperatures in all basins to ',num2str(defaultT), ' degree C.'));
    end
end
if ~isfield(PICO_opts,'Sbasins')
    defaultS = 33.8; % psu
    PICO_opts.Sbasins = zeros(Nbasins,1)+defaultS;
    if PICO_opts.InfoLevel>0
        warning(strcat('Ocean input vector is missing, setting salinities in all basins to ',num2str(defaultS), ' psu.'));
    end
end

if numel(PICO_opts.Sbasins)~=Nbasins
    error('Length of Salinity vector must equal number of basins');
end

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