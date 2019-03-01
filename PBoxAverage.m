function avg = PBoxAverage(MUA,Nind,Areas,Fin)

% Nind = PBOX==ID & ShelfID == SID;

pctBox = sum(Nind(MUA.connectivity),2)./3; % fraction of nodes in each element belonging to this box
BoxArea = sum(Areas.*pctBox);


msk = double(Nind(MUA.connectivity)); % Nele x Nnod where 1 indicates the node is within the box
msk(msk==0) = NaN; %mask out nodes that are in the element but not in the box so these don't count towards the box average
vals = nanmean(msk.*Fin(MUA.connectivity),2); % calculate the mean value in each element, not including nodes that are not within the box


avg = nansum(vals.*Areas.*pctBox)./BoxArea;
