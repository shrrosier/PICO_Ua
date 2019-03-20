%% Run Pico for the test example

%% Load data

load ResultsFiles/0000001-RegressionTest.mat;
%%
addpath(genpath('../')); % PATH TO PICO

%%
PrintInfoAboutElementsSizes(CtrlVar, MUA)

x = MUA.coordinates(:,1); 
y = MUA.coordinates(:,2);


%%
% Run PICO:
PICO_opts = struct;
PICO_opts.minArea = 1000;
PICO_opts.algorithm = 'watershed'; % PICO_opts.algorithm = 'polygon';
PICO_opts.C1 = 1e6;
PICO_opts.gamTstar = 2e-5;
PICO_opts.nmax = 5;
PICO_opts.SmallShelfMelt = 0;
PICO_opts.PICOres = 6000;
PICO_opts.Tbasins = [-1.76;-1.66;-1.65;-1.58;-1.51;-1.73;-1.68;-0.73;-1.61;-1.30;-1.83;-1.58;-0.36;0.47;1.04;1.17;0.23;-1.23;-1.80];
PICO_opts.Sbasins = [34.82;34.70;34.48;34.49;34.5;34.70;34.65;34.73;34.75;34.84;34.95;34.79;34.58;34.73;34.86;34.84;34.70;34.76;34.84];
PICO_opts.MeshBoundaryCoordinates = 'MeshBoundaryCoordinates.mat'

%tic

[Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(CtrlVar, MUA, GF, F.h, median(F.rho),F.rhow, PICO_opts);

%toc


save ResultsFiles/0000001-RegressionTest_PICO.mat PICO_opts Mk ShelfID T0 S0 Tkm Skm q PBOX Ak;

%%

do_plots = 0;

if do_plots
    %% Create figures:

    figure;
    PlotMeshScalarVariable(CtrlVar, MUA, PBOX)

    %% log melting
    decimals = 2; % number of decimal places that should be displayed +1

    Mk_log = sign(Mk).*log10(abs(Mk)*10^decimals);

    %
    figure; hold all;
    PlotMeshScalarVariable(CtrlVar, MUA, Mk_log)
    cbar = colorbar; caxis([-log10(0.1*10^decimals) log10(30*10^decimals)]);
    cbar.Label.String = 'sub-shelf melting (m/a)';
    cbar.TickLabels = {'-0.1' ,'0', '0.1', '1', '5', '10', '30'};
    cbar.Ticks = [-log10(0.1*10^decimals) 0  log10(0.1*10^decimals)  log10(1*10^decimals) log10(5*10^decimals) log10(10*10^decimals) log10(30*10^decimals)];
    PlotGroundingLines(CtrlVar, MUA, GF); PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar);

    %% aggregate melt rates over the ice shelves 

    ShelfID_per_ele = mean(ShelfID(MUA.connectivity),2);
    Int=FEintegrate2D([],MUA,Mk); % integarte melt rate over elements
    Areas = TriAreaFE(MUA.coordinates,MUA.connectivity); % get the area of each triangle

    average_melting_per_shelf = zeros(number_of_shelves,1);
    average_x_loc_per_shelf   = zeros(number_of_shelves,1);
    average_y_loc_per_shelf   = zeros(number_of_shelves,1);

    for shelf_i=1:number_of_shelves

        average_x_loc_per_shelf(shelf_i) = mean(x(ShelfID==shelf_i));
        average_y_loc_per_shelf(shelf_i) = mean(y(ShelfID==shelf_i));

        average_melting_per_shelf(shelf_i) = sum(Int(ShelfID_per_ele==shelf_i))/sum(Areas(ShelfID_per_ele==shelf_i));     

        % plot average melt rates 
        text(average_x_loc_per_shelf(shelf_i)/1000,average_y_loc_per_shelf(shelf_i)/1000, num2str(round(average_melting_per_shelf(shelf_i),2)) );

        % plot T0:
        %text(average_x_loc_per_shelf(shelf_i)/1000,average_y_loc_per_shelf(shelf_i)/1000, num2str(round(T0(shelf_i),2)) );

    end

    %%
    export_fig('Pico_melt_rates_polygon', '-dpdf', '-r400')

end;