
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
PICO_opts.algorithm = 'watershed';%'polygon';
PICO_opts.C1 = 1e6;
PICO_opts.gamTstar = 2e-5;
PICO_opts.nmax = 5;
PICO_opts.SmallShelfMelt = 0;
PICO_opts.PICOres = 2000;
PICO_opts.BasinsFile = 'BasinsInterpolant.mat';
PICO_opts.FloatingCriteria = 'GLthreshold'; %'GLthreshold' or 'StrictDownstream'

PICO_opts.Sbasins = [34.6505;34.5273;34.3222;34.3259;34.3297;34.5315;34.4819;34.5666;34.5766;34.6677;34.7822;34.6254;34.4107;34.5526;34.6902;34.6668;34.5339;34.5849;34.6644];
PICO_opts.Tbasins = [-1.75725;-1.65931;-1.58212;-1.54757;-1.51301;-1.72267;-1.69117;-0.67609;-1.61561;-1.30789;-1.83764;-1.5798;-0.368398;0.45505;1.04046;1.17196;0.229878;-1.23091;-1.79334];

PICO_opts.MeshBoundaryCoordinates = MeshBoundaryCoordinates;

tic

[Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(CtrlVarInRestartFile,MUA,GF,F.h,median(F.rho),F.rhow,PICO_opts);

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

ShelfID_per_ele = nanmean(ShelfID(MUA.connectivity),2);
Int=FEintegrate2D([],MUA,Mk); % integrate melt rate over elements
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

