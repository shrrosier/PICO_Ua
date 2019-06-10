function PICO_opts = PICO_DefaultParameters(MUA)
x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);

if ~exist('PICO_opts','var')
    warning('PICO_opts undefined, using default values... ARE YOU SURE YOU WANT TO DO THIS?');
    PICO_opts = struct;
end
if ~isfield(PICO_opts,'InfoLevel')
    fprintf('PICO InfoLevel undefined, setting PICO_opts.InfoLevel to maximum value\n');
    PICO_opts.InfoLevel = 100;
end
if ~isfield(PICO_opts,'algorithm')
    PICO_opts.algorithm = 'watershed';
end
if ~isfield(PICO_opts,'FloatingCriteria')
    PICO_opts.FloatingCriteria = 'GLthreshold';
end
if ~isfield(PICO_opts,'C1')
    PICO_opts.C1 = 1e6;
end
if ~isfield(PICO_opts,'gamTstar')
    PICO_opts.gamTstar = 2e-5;
end
if ~isfield(PICO_opts,'nmax')
    PICO_opts.nmax = 5;
end
if ~isfield(PICO_opts,'minArea')
    PICO_opts.minArea = 2e9;
end
if ~isfield(PICO_opts,'minNumShelf')
    PICO_opts.minNumShelf = 20;
end
if ~isfield(PICO_opts,'SmallShelfMelt')
    PICO_opts.SmallShelfMelt = 0;
    if PICO_opts.InfoLevel>0
        warning('Applying default zero melt to floating nodes outside of delineated shelves, change this in PICO_opts.SmallShelfMelt');
    end
end
if ~isfield(PICO_opts,'BasinsFile')
    Fbasins = scatteredInterpolant(x,y,ones(MUA.Nnodes,1));
    save('DefaultBasinsInterpolant.mat', 'Fbasins');
    PICO_opts.BasinsFile = 'DefaultBasinsInterpolant.mat';
    PICO_opts.Nbasins = 1;
    if PICO_opts.InfoLevel>0
        warning('Basins file missing, setting everything to one basin');
    end
end
if ~isfield(PICO_opts, 'Nbasins')
    load(PICO_opts.BasinsFile)
    PICO_opts.Nbasins = max(Fbasins(x,y));
end
if ~isfield(PICO_opts,'Tbasins')
    defaultT = -1.8; % deg C
    PICO_opts.Tbasins = zeros(PICO_opts.Nbasins,1)+defaultT;
    if PICO_opts.InfoLevel>0
        warning(strcat('Ocean input vector is missing, setting temperatures in all basins to ',num2str(defaultT), ' degree C.'));
    end
end
if ~isfield(PICO_opts,'Sbasins')
    defaultS = 33.8; % psu
    PICO_opts.Sbasins = zeros(PICO_opts.Nbasins,1)+defaultS;
    if PICO_opts.InfoLevel>0
        warning(strcat('Ocean input vector is missing, setting salinities in all basins to ',num2str(defaultS), ' psu.'));
    end
end
if numel(rhoi) > 1
    if PICO_opts.InfoLevel>0
        warning('non scalar rhoi, averaging...');
    end
    rhoi = mean(rhoi);
end
if numel(rhow) > 1
    if PICO_opts.InfoLevel>0
        warning('non scalar rhow, averaging...');
    end
    rhow = mean(rhow);
end
switch PICO_opts.algorithm
    case 'watershed'
        if ~isfield(PICO_opts,'PICOres')
            % try to make a sensible choice about resolution vs speed
            if MUA.Nnodes > 200e3
                PICO_opts.PICOres = 6000;
            elseif MUA.Nnodes > 100e3
                PICO_opts.PICOres = 3000;
            else
                PICO_opts.PICOres = 1000;
            end
            if PICO_opts.InfoLevel>0
                warning('Using default resolution, change this in PICO_opts.PICOres');
            end
        end
        if PICO_opts.InfoLevel>1
            fprintf('Using watershed algorithm to deliniate ice shelves...\n');
        end
        
        [ShelfID,PBOX,Ak,floating] = IdentifyIceShelvesWatershedOption(CtrlVar,MUA,GF,PICO_opts.PICOres,PICO_opts.minArea,PICO_opts.minNumShelf,PICO_opts.nmax,PICO_opts.FloatingCriteria);
        
    case 'polygon'
        if ~isfield(PICO_opts,'MeshBoundaryCoordinates')
            error('PICO_opts.MeshBoundaryCoordinates must be defined for the polygon option');
        end
        if ~isfield(PICO_opts,'persistentBC')
            PICO_opts.persistentBC = 0;
        end
        if PICO_opts.InfoLevel>1
            fprintf('Using polygon algorithm to deliniate ice shelves...\n');
        end
        [ShelfID,PBOX,Ak,floating] = IdentifyIceShelvesPolygonOption(CtrlVar,MUA,GF,PICO_opts);
    otherwise
        error('Invalid algorithm, choose either "watershed" or "polygon"');
end




end