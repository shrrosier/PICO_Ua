function [Mk,ShelfID,T0,S0,Tkm,Skm,q] = PICO_driver(CtrlVar,MUA,GF,b,rhoi,PICO_opts)
%
% PICO melt rate parameterisation v0.4
% Outputs melt rate in units of metres per year
% NOTE:
% THE WATERSHED OPTION MAY BE FASTER BUT REQUIRES THE IMAGE PROCESSING TOOLBOX
%
% Exampple Usage:
%   [Mk,ShelfID,T0,S0] = PICO_driver(CtrlVar,MUA,GF,b,rho,PICO_opts)
%
% Inputs:
%
%   CtrlVar: Control Variables (Ua structure)
%   MUA: Mesh (Ua structure)
%   GF: Grounded/floating mask (Ua structure)
%   b: ice base (MUA.Nnodes x 1)
%   rho: ice density (MUA.Nnodes x 1)
%   PICO_opts: structure containing various PICO options, listed below:
%
%   - algorithm: either 'watershed' or 'polygon'
%   - PICOres: Resolution to use when searching for ice shelves i.e. 10000
%   (only used in watershed algorithm)
%   - SmallShelfMelt: melt to floating regions not included in PICO boxes i.e. 0
%   - nmax: maximum number of PICO boxes (optional, default = 5)
%   - minArea: minimum area for an ice shelf to be defined in PICO
%   - minNumShelf: minimum number of floating nodes for an ice shelf to be
%   defined in PICO
%   - C: Overturning strength (optional, default = 1e6)
%   - gamTstar: turbulent temp. exch. coeff. (optional, default = 2e-5)
%   - BasinsFile: .mat file containing scattered interpolant of Basin IDs
%   for you domain
%   - Tbasins: for each basin ID, the corresponding temperature T0
%   - Sbasins: for each basin ID, the corresponding salinity S0
%   - MeshBoundaryCoordinates: Nnodes x 2 matrix containing the mesh boundary
%   coordinates, as defined in Ua2DInitialUserInput (for polygon option)
%
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);

if ~exist('PICO_opts','var')
    PICO_opts = struct;
end
if ~isfield(PICO_opts,'algorithm')
    PICO_opts.algorithm = 'watershed';
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
    warning('Applying default zero melt to floating nodes outside of floating, change this in PICO_opts.SmallShelfMelt');
end
if ~isfield(PICO_opts,'BasinsFile')
    Fbasins = scatteredInterpolant(x,y,ones(MUA.Nnodes,1));
    save('DefaultBasinsInterpolant.mat', 'Fbasins');
    PICO_opts.BasinsFile = 'DefaultBasinsInterpolant.mat';
    PICO_opts.Nbasins = 1;
    warning('Basins file missing, setting everything to one basin');
end
if ~isfield(PICO_opts, 'Nbasins')
    load(PICO_opts.BasinsFile)
    PICO_opts.Nbasins = max(Fbasins(x,y));
    fprintf(strcat('Using '));
end
if ~isfield(PICO_opts,'Tbasins')
    defaultT = -1.8; % deg C
    PICO_opts.Tbasins = zeros(PICO_opts.Nbasins,1)+defaultT;
    warning(strcat('Ocean input vector is missing, setting temperatures in all basins to ',num2str(defaultT), ' degree C.'));
end
if ~isfield(PICO_opts,'Sbasins')
    defaultS = 33.8; % psu
    PICO_opts.Sbasins = zeros(PICO_opts.Nbasins,1)+defaultS;
    warning(strcat('Ocean input vector is missing, setting salinities in all basins to ',num2str(defaultS), ' psu.'));
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
            warning('Using default resolution, change this in PICO_opts.PICOres');
        end
        
        [ShelfID,PBOX,Ak,floating] = IdentifyIceShelvesWatershedOption(CtrlVar,MUA,GF,PICO_opts.PICOres,PICO_opts.minArea,PICO_opts.minNumShelf,PICO_opts.nmax);
    case 'polygon'
        if ~isfield(PICO_opts,'MeshBoundaryCoordinates')
            error('PICO_opts.MeshBoundaryCoordinates must be defined for the polygon option');
        end
        [ShelfID,PBOX,Ak,floating] = IdentifyIceShelvesPolygonOption(CtrlVar,MUA,GF,PICO_opts.minArea,PICO_opts.minNumShelf,PICO_opts.nmax,PICO_opts.MeshBoundaryCoordinates);
    otherwise
        error('Invalid algorithm, choose either "watershed" or "polygon"');
end

% ========================= input from box b0 =====================

[T0,S0] = GetOceanInputFromBox0(CtrlVar,MUA,GF,ShelfID,PICO_opts);

% ========================= physics ===============================
% options to specify for PICO: C, gamTstar

a1 = -0.0572; % deg C PSU-1
b1 = 0.0788; % deg C
c1 = 7.77e-8; % deg C Pa-1
alph = 7.5e-5; % deg C-1
beta = 7.7e-4; % PSU-1
rhostar = 1033;
L = 3.34e5; %J Kg-1
cp = 3974; %J Kg-1 degC-1

C1 = PICO_opts.C1;
gamTstar = PICO_opts.gamTstar;
nmax = PICO_opts.nmax;

% fix this... only need them for mu, calculate over each box?
rhow = 1030;
rho = 910;

Sk = zeros(MUA.Nnodes,1);
Tk = zeros(MUA.Nnodes,1);
Tkm = zeros(max(ShelfID),nmax);
Skm = zeros(max(ShelfID),nmax);
Tstar = zeros(MUA.Nnodes,1);


mu = rho./rhow;
lambda = L/cp;
s1 = S0./(mu*lambda);
gk = Ak.*gamTstar;
pk = -b.*9.81.*rhoi;
gk2 = gk./(mu*lambda);


% ========================= box 1 ===============================

pcoeff = get_p(gk(:,1),s1);

for ii = 1:max(ShelfID)
    
    ind = ShelfID==ii & PBOX == 1;
    
    Tstar(ind) = calc_tstar(S0(ii),T0(ii),pk(ind));
    if any(Tstar > 0) % equivalent to: if T0 < T_pmp
        T_pmp = Tstar+T0(ii);
         % ensure that temperature input to box 1 at least at pressure melting point
         T0(ii) = max(T_pmp(ind)+0.001,T0(ii));
        Tstar(ind) = calc_tstar(S0(ii),T0(ii),pk(ind));
    end
    q1 = pcoeff(ii).*Tstar;
    
    DD = .25.*pcoeff(ii).^2 - q1;
    if any(DD<0)
        error('something went wrong here - ask Ronja!');
    end
    
    Tk(ind) = T0(ii) - (-.5.*pcoeff(ii) + sqrt(DD(ind)));
    Sk(ind) = S0(ii) - (S0(ii)/(mu*lambda)) * (T0(ii) - Tk(ind));
    
    Skm(ii,1) = mean(Sk(ind));
    Tkm(ii,1) = mean(Tk(ind));
    
end

q = C1*rhostar*(beta*(S0-Skm(:,1)) - alph*(T0-Tkm(:,1)));

% ========================= box k ===============================

for ii = 1:max(ShelfID)
    for jj = 2:nmax
        
        ind = ShelfID==ii & PBOX == jj;        
        
        
        Tstar = calc_tstar(Skm(ii,jj-1),Tkm(ii,jj-1),pk(ind)); % calculate with Sk and Tk of box-1 but using local pk
        
        Tk(ind) = Tkm(ii,jj-1) + (gk(ii,jj).*Tstar)./(q(ii) + gk(ii,jj) - gk2(ii,jj).*a1.*Skm(ii,jj-1));
        Tkm(ii,jj) = mean(Tk(ind));
        Sk(ind) = Skm(ii,jj-1) - Skm(ii,jj-1).*(Tkm(ii,jj-1)-Tkm(ii,jj))./(mu*lambda);
        Skm(ii,jj) = mean(Sk(ind));
    end
end

Mk_ms = (-gamTstar./(mu*lambda)).*(a1.*Sk + b1 - c1.*pk - Tk);
Mk = Mk_ms .* 86400 .* 365.25;
Mk(~floating) = 0;
Mk(isnan(ShelfID) & floating) = PICO_opts.SmallShelfMelt;

    function p_coeff = get_p(g1,s1)
        p_coeff = g1./(C1*rhostar*(beta*s1 - alph));
    end

    function Tpm = calc_pmp(S,p)
        Tpm = a1.*S + b1 - c1.*p;
    end

    function Tstar = calc_tstar(S,T,p)
        tpm = calc_pmp(S,p);
        Tstar = tpm - T;
    end

end
