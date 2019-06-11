function [Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(UserVar,CtrlVar,MUA,GF,h,rhoi,rhow,PICO_opts)
%
% Usage:
% [Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(CtrlVar,MUA,GF,h,rhoi,rhow,PICO_opts)
%
% Mk is melt rate in metres per year
%
% optional outputs:
%
% ShelfID: Shelf ID number of each node
% T0,S0: ambient basin T and S of each node
% Tkm, Skm: calculated T and S of each node
% q: overturning
% PBOX: Box ID number of each node
% Ak: Areas of each Shelf
%
% type 'help PICO' for more details
%

PICO_opts = PICO_DefaultParameters(MUA,PICO_opts);

if numel(rhoi) > 0
    if PICO_opts.InfoLevel>0
        warning('non scalar rhoi, averaging...');
    end
    rhoi = mean(rhoi);
end
if numel(rhow) > 0
    if PICO_opts.InfoLevel>0
        warning('non scalar rhow, averaging...');
    end
    rhow = mean(rhow);
end

switch PICO_opts.algorithm
    case 'watershed'
        
        if PICO_opts.InfoLevel>1
            fprintf('Using watershed algorithm to deliniate ice shelves...\n');
        end
        
        [ShelfID,PBOX,Ak,floating] = PICO_IdentifyIceShelvesWatershedOption(UserVar,CtrlVar,MUA,GF,PICO_opts.PICOres,PICO_opts.minArea,PICO_opts.minNumShelf,PICO_opts.nmax,PICO_opts.FloatingCriteria);
        
    case 'polygon'
        if PICO_opts.InfoLevel>1
            fprintf('Using polygon algorithm to deliniate ice shelves...\n');
        end
        [ShelfID,PBOX,Ak,floating] = PICO_IdentifyIceShelvesPolygonOption(UserVar,CtrlVar,MUA,GF,PICO_opts);
        
    case 'oneshelf'
        if PICO_opts.InfoLevel>1
            fprintf('using oneshelf option, ONLY DO THIS IF YOURE CERTAIN YOU WILL ONLY EVER HAVE ONE ICE SHELF...\n');
        end
        [ShelfID,PBOX,Ak,floating] = PICO_MakeBoxesOneShelf(UserVar,CtrlVar,MUA,GF,PICO_opts);
        
    otherwise
        error('Invalid algorithm, choose either "watershed" or "polygon"');
end

if PICO_opts.InfoLevel>0
    fprintf('%s algorithm completed, found %2i ice shelves in the domain\n',PICO_opts.algorithm,max(ShelfID));
end


if PICO_opts.InfoLevel>10
    
    fprintf('Plotting delineated ice shelves...\n');
    x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
    figure(809);
    PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k');
    hold on;
    PlotMeshScalarVariable(CtrlVar, MUA, ShelfID);
    colormap('lines')
    colorbar off
    
    for shelf_i=1:max(ShelfID)
        average_x_loc_per_shelf(shelf_i) = mean(x(ShelfID==shelf_i));
        average_y_loc_per_shelf(shelf_i) = mean(y(ShelfID==shelf_i));
        text(mean(x(ShelfID==shelf_i))/1000,mean(y(ShelfID==shelf_i))/1000,num2str(shelf_i));
    end
    title('Shelf IDs');
    
    fprintf('Plotting ice shelf boxes...\n');
    figure(811);
    PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k');
    hold on;
    tempBoxID = PBOX; tempBoxID(tempBoxID==0) = NaN;
    PlotMeshScalarVariable(CtrlVar, MUA, tempBoxID);
    colormap(summer(PICO_opts.nmax));
    title('Box Numbers');
end

% ========================= input from box b0 =====================
if PICO_opts.InfoLevel>1
    fprintf('Calculating temperature and salinity for box 0...\n');
end
[T0,S0] = PICO_GetOceanInputFromBox0(UserVar,CtrlVar,MUA,GF,ShelfID,PICO_opts);

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

Sk = zeros(MUA.Nnodes,1);
Tk = zeros(MUA.Nnodes,1);
Tkm = zeros(max(ShelfID),nmax);
Skm = zeros(max(ShelfID),nmax);
Tstar = zeros(MUA.Nnodes,1);

mu = rhoi./rhow;
lambda = L/cp;
s1 = S0./(mu*lambda);
gk = Ak.*gamTstar;
pk = h.*9.81.*rhoi;
gk2 = gk./(mu*lambda);


% ========================= box 1 ===============================

if PICO_opts.InfoLevel>1
    fprintf('Calculating temperature and salinity for box 1...\n');
end

Areas=TriAreaFE(MUA.coordinates,MUA.connectivity);
ShelfIDEle = nanmean(ShelfID(MUA.connectivity),2);

% ------------------------------
% pctBox and mskBox are needed
% for averaging values in boxes
for ii = 1:nmax
    Nind = PBOX ==ii;
    msk = double(Nind(MUA.connectivity));
    pctBox{ii} = sum(msk,2)./3;
    msk(msk==0) = nan;
    mskBox{ii} = msk;
end
% ------------------------------

pcoeff = get_p(gk(:,1),s1);

for ii = 1:max(ShelfID)
    
    ind = ShelfID==ii & PBOX == 1;
    indE = ShelfIDEle==ii;
    
    pb = pctBox{1}(indE);
    mb = mskBox{1}(indE,:);
    ar = Areas(indE);
    BA = sum(pb.*ar);
    
    
    Tstar(ind) = calc_tstar(S0(ii),T0(ii),pk(ind));
    if any(Tstar(ind) > 0) % equivalent to: if T0 < T_pmp
        % ensure that temperature input to box 1 is at least at pressure melting point
        J = Tstar > 0 & ind;
        T_pmp = Tstar+T0(ii);
        Tstar(J) = calc_tstar(S0(ii),T_pmp(J)+0.001,pk(J));
    end
    q1 = pcoeff(ii).*Tstar;
    
    DD = .25.*pcoeff(ii).^2 - q1;
    if any(DD<0)
        error('something went wrong here!');
    end
    
    Tk(ind) = T0(ii) - (-.5.*pcoeff(ii) + sqrt(DD(ind)));
    Sk(ind) = S0(ii) - (S0(ii)/(mu*lambda)) * (T0(ii) - Tk(ind));
    
    vals = nanmean(mskBox{1}.*Sk(MUA.connectivity),2);
    Skm(ii,1) = nansum(vals(indE).*ar.*pb)./BA;
    
    vals = nanmean(mskBox{1}.*Tk(MUA.connectivity),2);
    Tkm(ii,1) = nansum(vals(indE).*ar.*pb)./BA;
    
end

q = C1*rhostar*(beta*(S0-Skm(:,1)) - alph*(T0-Tkm(:,1)));

% ========================= box k ===============================

if PICO_opts.InfoLevel>10
    fprintf('Calculating temperature and salinity for ice shelf ');
end

for ii = 1:max(ShelfID)
    
    indE = ShelfIDEle==ii;
    ar = Areas(indE);
    
    for jj = 2:nmax
        
        if Ak(ii,jj) == 0
            continue
        end
        
        ind = ShelfID==ii & PBOX == jj;
        
        pb = pctBox{jj}(indE);
        mb = mskBox{jj}(indE,:);
        BA = sum(pb.*ar);
        
        Tstar = calc_tstar(Skm(ii,jj-1),Tkm(ii,jj-1),pk(ind)); % calculate with Sk and Tk of box-1 but using local pk
        
        Tk(ind) = Tkm(ii,jj-1) + (gk(ii,jj).*Tstar)./(q(ii) + gk(ii,jj) - gk2(ii,jj).*a1.*Skm(ii,jj-1));
        vals = nanmean(mskBox{jj}.*Tk(MUA.connectivity),2);
        Tkm(ii,jj) = nansum(vals(indE).*ar.*pb)./BA;
        
        Sk(ind) = Skm(ii,jj-1) - Skm(ii,jj-1).*(Tkm(ii,jj-1)-Tkm(ii,jj))./(mu*lambda);
        vals = nanmean(mskBox{jj}.*Sk(MUA.connectivity),2);
        Skm(ii,jj) = nansum(vals(indE).*ar.*pb)./BA;
        
        
    end
    
    if PICO_opts.InfoLevel>10
        fprintf(strcat(num2str(ii),', '))
    end
    
end
% calculate melt rate (Mk)

if PICO_opts.InfoLevel>1
    fprintf('\nCalculating basal melt rates...\n');
end

Mk_ms = (-gamTstar./(mu*lambda)).*(a1.*Sk + b1 - c1.*pk - Tk); % m per s
Mk = Mk_ms .* 86400 .* 365.25; % m per a
Mk(~floating) = 0;
Mk(isnan(ShelfID) & floating) = PICO_opts.SmallShelfMelt;

if PICO_opts.InfoLevel>10
    fprintf('\n');
    fprintf('----------|----------------|--------|---------|--------|----------|----------|----------|\n');
    fprintf('Shelf no. |   Shelf Area   |  T0    |   S0    | Tfront |  Sfront  |  max ab  |  min ab  |\n');
    fprintf('----------| ---------------|--------|---------|--------|----------|----------|----------|\n');
    for ii = 1:max(ShelfID)
        fprintf('    %2i    | %14.7g | %6.3g | %7.4g | %6.3g | %7.4g  | %8.5g | %8.5g |\n',ii,sum(Ak(ii,:)),T0(ii),S0(ii),Tkm(max(PBOX(ShelfID==ii)))...
            ,Skm(max(PBOX(ShelfID==ii))),max(Mk(ShelfID==ii)),min(Mk(ShelfID==ii)));
    end
    fprintf('----------|----------------|--------|---------|--------|----------|----------|----------|\n');
elseif PICO_opts.InfoLevel>0
    fprintf('PICO run complete, ice shelves have max, mean and min melt rates of %-g,  %-g,  %-g, respectively. \n',max(Mk),mean(Mk),min(Mk));
end

if PICO_opts.InfoLevel>10
    fprintf('Plotting PICO melt rates for all ice shelves...\n');
    
    decimals = 2; % number of decimal places that should be displayed +1
    Mk_log = sign(Mk).*log10(abs(Mk)*10^decimals);
    figure(821); hold all;
    PlotMeshScalarVariable(CtrlVar, MUA, Mk_log);
    cbar = colorbar; caxis([-log10(0.1*10^decimals) log10(30*10^decimals)]);
    cbar.Label.String = 'sub-shelf melting (m/a)';
    cbar.TickLabels = {'-0.1' ,'0', '0.1', '1', '5', '10', '30'};
    cbar.Ticks = [-log10(0.1*10^decimals) 0  log10(0.1*10^decimals)  log10(1*10^decimals) log10(5*10^decimals) log10(10*10^decimals) log10(30*10^decimals)];
    PlotGroundingLines(CtrlVar, MUA, GF); PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar);
    
    % aggregate melt rates over the ice shelves (first guess, ignore the area of the elements)
    
    ShelfID_per_ele = nanmean(ShelfID(MUA.connectivity),2);
    Int=FEintegrate2D([],MUA,Mk); % integrate melt rate over elements
    Areas = TriAreaFE(MUA.coordinates,MUA.connectivity); % get the area of each triangle
    number_of_shelves = max(ShelfID);
    
    average_melting_per_shelf = zeros(number_of_shelves,1);
    average_x_loc_per_shelf   = zeros(number_of_shelves,1);
    average_y_loc_per_shelf   = zeros(number_of_shelves,1);
    
    for shelf_i=1:number_of_shelves
        average_x_loc_per_shelf(shelf_i) = mean(x(ShelfID==shelf_i));
        average_y_loc_per_shelf(shelf_i) = mean(y(ShelfID==shelf_i));
        average_melting_per_shelf(shelf_i) = sum(Int(ShelfID_per_ele==shelf_i))/sum(Areas(ShelfID_per_ele==shelf_i));
        text(average_x_loc_per_shelf(shelf_i)/1000,average_y_loc_per_shelf(shelf_i)/1000, num2str(round(average_melting_per_shelf(shelf_i),2)) );
    end
    
    title('PICO melt rates');
end


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
