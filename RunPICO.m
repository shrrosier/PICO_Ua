
%% Load data

load Ua2D_Restartfile.mat;
load MeshBoundaryCoordinates.mat

%%
PrintInfoAboutElementsSizes(CtrlVarInRestartFile, MUA)

x = MUA.coordinates(:,1); 
y = MUA.coordinates(:,2);


%%
% Run PICO:
PICO_opts = struct;
PICO_opts.algorithm = 'polygon';%'polygon'; %'watershed';
PICO_opts.C1 = 1e6;
PICO_opts.gamTstar = 2e-5;
PICO_opts.nmax = 5;
PICO_opts.SmallShelfMelt = 0;
PICO_opts.PICOres = 6000;
PICO_opts.BasinsFile = 'BasinsInterpolant.mat';
PICO_opts.FloatingCriteria = 'GLthreshold'; %'GLthreshold' or 'StrictDownstream'

PICO_opts.Sbasins = [34.6505; 34.5273; 34.3222; 34.3259; 34.3297; 34.5315; 34.4819; 34.5666; 34.5766; 34.6677; 34.7822; 34.6254; 34.4107; 34.5526; 34.6902; 34.6668; 34.5339; 34.5849; 34.6644];
PICO_opts.Tbasins = [-1.76;-1.66;-1.65;-1.58;-1.51;-1.73;-1.68;-0.73;-1.61;-1.30;-1.83;-1.58;-0.36;0.47;1.04;1.17;0.23;-1.23;-1.80];
% PICO_opts.Sbasins = [34.82;34.70;34.48;34.49;34.5;34.70;34.65;34.73;34.75;34.84;34.95;34.79;34.58;34.73;34.86;34.84;34.70;34.76;34.84];
PICO_opts.MeshBoundaryCoordinates = MeshBoundaryCoordinates;

tic

[Mk,ShelfID,T0,S0,Tk,Sk,q] = PICO_driver(CtrlVarInRestartFile,MUA,GF,F.h,F.rho, PICO_opts);

toc


% Infos from profiling: 
% Running "watershed", res 1000: Elapsed time is 45.102403 seconds.
% Running "watershed", res 3000: Elapsed time is  6.471619 seconds.
% Running "watershed", res 6000: Elapsed time is  2.838979 seconds.
% Running 


%% Create figures:

decimals = 2; % number of decimal places that should be displayed +1

Mk_log = sign(Mk).*log10(abs(Mk)*10^decimals);

%
figure; hold all;
PlotMeshScalarVariable(CtrlVarInRestartFile, MUA, Mk_log)
cbar = colorbar; caxis([-log10(0.1*10^decimals) log10(30*10^decimals)]);
cbar.Label.String = 'sub-shelf melting (m/a)';
cbar.TickLabels = {'-0.1' ,'0', '0.1', '1', '5', '10', '30'};
cbar.Ticks = [-log10(0.1*10^decimals) 0  log10(0.1*10^decimals)  log10(1*10^decimals) log10(5*10^decimals) log10(10*10^decimals) log10(30*10^decimals)];
PlotGroundingLines(CtrlVarInRestartFile, MUA, GF); PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVarInRestartFile);

%% aggregate melt rates over the ice shelves (first guess, ignore the area of the elements)

ShelfID_per_ele = mean(ShelfID(MUA.connectivity),2);
Int=FEintegrate2D([],MUA,Mk); % integarte melt rate over elements
Areas = TriAreaFE(MUA.coordinates,MUA.connectivity); % get the area of each triangle

% FIXME: not aggregated per region!!
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
    
end

% %%
% export_fig('Pico_melt_rates_polygon', '-dpdf', '-r400')

