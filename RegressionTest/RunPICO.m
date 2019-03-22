function Passed=RunPICO(option)


load ./ResultsFiles/0000001-RegressionTest.mat;
load MeshBoundaryCoordinates.mat


PrintInfoAboutElementsSizes(CtrlVar, MUA)


%%
% Run PICO:

PICO_opts = struct;
PICO_opts.minArea = 1000; % cut-off area for ice shelves
PICO_opts.algorithm = option; % watershed or polygon
PICO_opts.C1 = 1e6; % overturning coefficient
PICO_opts.gamTstar = 2e-5; % heat exchange velocity
PICO_opts.nmax = 5; % maximum number of boxes (used for large ice shelves)
PICO_opts.SmallShelfMelt = 0; % melting in shelves not dealt with in PICO
PICO_opts.PICOres = 6000;
PICO_opts.Tbasins = -1.5; % input temperature in box 0
PICO_opts.Sbasins = 34.5; % input salinity in box 0
PICO_opts.MeshBoundaryCoordinates = MeshBoundaryCoordinates;

%tic

[Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(CtrlVar, MUA, GF, F.h, median(F.rho),F.rhow, PICO_opts);

%toc

save ResultsFiles/0000001-RegressionTest_PICO.mat PICO_opts Mk ShelfID T0 S0 Tkm Skm q PBOX Ak;

%%

do_plots = 1;

if do_plots
    

    x = MUA.coordinates(:,1); 
    y = MUA.coordinates(:,2);
    
    number_of_shelves = PICO_opts.nmax;
    
    % boxes
    figure;
    PlotMeshScalarVariable(CtrlVar, MUA, PBOX)
    export_fig(strcat('Pico_boxes_',PICO_opts.algorithm), '-dpdf', '-r400')
    
    % log melting
    decimals = 2; % number of decimal places that should be displayed +1

    Mk_log = sign(Mk).*log10(abs(Mk)*10^decimals);

    figure; hold all;
    PlotMeshScalarVariable(CtrlVar, MUA, Mk_log)
    cbar = colorbar; caxis([-log10(0.1*10^decimals) log10(30*10^decimals)]);
    cbar.Label.String = 'sub-shelf melting (m/a)';
    cbar.TickLabels = {'-0.1' ,'0', '0.1', '1', '5', '10', '30'};
    cbar.Ticks = [-log10(0.1*10^decimals) 0  log10(0.1*10^decimals)  log10(1*10^decimals) log10(5*10^decimals) log10(10*10^decimals) log10(30*10^decimals)];
    PlotGroundingLines(CtrlVar, MUA, GF); PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar);

    %%
    export_fig(strcat('Pico_melt_rates_',PICO_opts.algorithm), '-dpdf', '-r400')

end;

end