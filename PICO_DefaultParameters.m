function PICO_opts = PICO_DefaultParameters(MUA,PICO_opts)
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
    if PICO_opts.InfoLevel>1
        fprintf('PICO algorithm undefined, defaulting to watershed\n');
    end
end
if ~isfield(PICO_opts,'FloatingCriteria')
    PICO_opts.FloatingCriteria = 'GLthreshold';
    if PICO_opts.InfoLevel>10
        fprintf('Floating Criterion undefined, using GLthreshold\n');
    end
end
if ~isfield(PICO_opts,'C1')
    PICO_opts.C1 = 1e6;
    if PICO_opts.InfoLevel>1
        warning('Overturning strength undefined, using C = 1e6\n');
    end
end
if ~isfield(PICO_opts,'gamTstar')
    PICO_opts.gamTstar = 2e-5;
    if PICO_opts.InfoLevel>1
        warning('Turbulent temperature exchange velocity undefined , using gamTstar = 2e-5\n');
    end
end
if ~isfield(PICO_opts,'nmax')
    PICO_opts.nmax = 5;
        if PICO_opts.InfoLevel>10
        fprintf('Using Default number of boxes (5)\n');
    end
end
if ~isfield(PICO_opts,'minArea')
    PICO_opts.minArea = 2e9;
    if PICO_opts.InfoLevel>1
        warning('Using default shelf area cutoff of 2e9, shelves smaller will be ignored!\n');
    end
end
if ~isfield(PICO_opts,'minNumShelf')
    PICO_opts.minNumShelf = 20;
        if PICO_opts.InfoLevel>1
        warning('Using default shelf size cutoff of 20, shelves smaller will be ignored!\n');
    end
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
    if PICO_opts.InfoLevel>0
        warning('Basins file missing, setting everything to one basin');
    end
end

if ~isfield(PICO_opts,'PICOres') && strcmp(PICO_opts.algorithm,'watershed')
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

if ~isfield(PICO_opts,'MeshBoundaryCoordinates')  && strcmp(PICO_opts.algorithm,'polygon')
    error('PICO_opts.MeshBoundaryCoordinates must be defined for the polygon option');
end

if ~isfield(PICO_opts,'persistentBC') && strcmp(PICO_opts.algorithm,'polygon')
    PICO_opts.persistentBC = 0;
end




end