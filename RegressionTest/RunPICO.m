%%
% Watershed: only shelf is identified as bad shelf
% Polygon: problem with MeshBoundaryCoordinates


%% Load data
clearvars
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
PICO_opts.algorithm = 'polygon';
PICO_opts.C1 = 1e6;
PICO_opts.gamTstar = 2e-5;
PICO_opts.nmax = 5;
PICO_opts.SmallShelfMelt = 0;
PICO_opts.PICOres = 500;
PICO_opts.Tbasins = -1.0;
PICO_opts.Sbasins = 34.5; 
load('MeshBoundaryCoordinates.mat');
PICO_opts.MeshBoundaryCoordinates = MeshBoundaryCoordinates;

%tic

[Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX] = PICO_driver(CtrlVar, MUA, GF, F.h, median(F.rho), PICO_opts);

%toc

% Infos from profiling: 
% Running "watershed", res 1000: 
% Running "watershed", res 3000: 
% Running "watershed", res 6000: 
% Running 
%%

Tkm
Skm
q
median(Mk(PBOX==1))
%% Create figures:

figure;
PlotMeshScalarVariable(CtrlVar, MUA, PBOX)


%% Old:

%% Calculate area of boxes -> are they equal?


%% 




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

%% aggregate melt rates over the ice shelves (first guess, ignore the area of the elements)

ShelfID_per_ele = mean(ShelfID(MUA.connectivity),2);
Int=FEintegrate2D([],MUA,Mk); % integarte melt rate over elements
Areas = TriAreaFE(MUA.coordinates,MUA.connectivity); % get the area of each triangle


number_of_shelves = max(ShelfID);

average_melting_per_shelf = zeros(number_of_shelves,1);
average_x_loc_per_shelf   = zeros(number_of_shelves,1);
average_y_loc_per_shelf   = zeros(number_of_shelves,1);

for shelf_i=1:number_of_shelves

    average_x_loc_per_shelf(shelf_i) = mean(x(ShelfID==shelf_i));
    average_y_loc_per_shelf(shelf_i) = mean(y(ShelfID==shelf_i));
    
    average_melting_per_shelf(shelf_i) = sum(Int(ShelfID_per_ele==shelf_i))/sum(Areas(ShelfID_per_ele==shelf_i));     
    
    % plot avreage melt rates 
    text(average_x_loc_per_shelf(shelf_i)/1000,average_y_loc_per_shelf(shelf_i)/1000, num2str(round(average_melting_per_shelf(shelf_i),2)) );
    
    % plot also T0:
    %text(average_x_loc_per_shelf(shelf_i)/1000,average_y_loc_per_shelf(shelf_i)/1000, num2str(round(T0(shelf_i),2)) );
    
end

%%
export_fig('Pico_melt_rates_polygon', '-dpdf', '-r400')

