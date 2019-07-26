function Ak = PICO_calc_areas(MUA, BoxID, ShelfID, nmax)

% this function takes connectivity, box # and shelf # to produce a matrix
% of areas with dimensions # shelves x # boxes
% The issue here is that many elements will contain nodes in different
% boxes or nodes that are both grounded and floating

[Areas,~,~,~]=TriAreaFE(MUA.coordinates,MUA.connectivity);

Areapct = cell(1,nmax); % the percent of every element that is contained in a box number given by the cell number ie cell one = box one
for ii = 1:nmax
    Nind = BoxID ==ii;
    msk = double(Nind(MUA.connectivity));
    pctBox = sum(msk,2)./3;
    Areapct{ii} = pctBox.*Areas;
end

ShelfIDEle = round(nanmean(double(ShelfID(MUA.connectivity)),2));
Ak = zeros(max(ShelfID),nmax);

for ii = 1:max(ShelfID)
    ind = ShelfIDEle==ii;
    for k = 1:nmax
        Ak(ii,k) = sum(Areapct{k}(ind));
    end
end

