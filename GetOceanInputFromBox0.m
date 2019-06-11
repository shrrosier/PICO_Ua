function [T0,S0]=GetOceanInputFromBox0(UserVar,CtrlVar,MUA,GF,ShelfID,PICO_opts)
% Calculate the ocean input for each ice shelf from the ocean input 
% values per basin, weighted with the area of the ice shelf within that
% basin

T0 = zeros(max(ShelfID),1);
S0 = zeros(max(ShelfID),1);

load(PICO_opts.BasinsFile); % Fbasins

basins = Fbasins(MUA.coordinates(:,1), MUA.coordinates(:,2));

MUA2 = MUA;
MUA2.nip = 1;
TriArea = FEintegrate2D([],MUA2,ones(MUA.Nnodes,1));

for ii = 1:max(ShelfID)
    Ak = zeros(PICO_opts.Nbasins,1);
    % calc the area of the ice shelf lying in each basin...
    for basin_i = 1:PICO_opts.Nbasins
        ind = ShelfID==ii & basins==basin_i;
        Ak(basin_i)=sum(TriArea(ind));
    end
    T0(ii) = dot(PICO_opts.Tbasins, 1/sum(Ak)*Ak);
    S0(ii) = dot(PICO_opts.Sbasins, 1/sum(Ak)*Ak);
end

end