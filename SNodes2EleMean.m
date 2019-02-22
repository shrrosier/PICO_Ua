function [MeanEleValues]=SNodes2EleMean(connectivity,NodalValues)
    
    %
    % [MeanEleValues]=Nodes2EleMean(connectivity,NodalValues)
    % Takes nodal values and returns average of nodal values for each element
    %
    [Nele,nod]=size(connectivity);
    MeanEleValues=max(reshape(NodalValues(connectivity),Nele,nod),[],2);
    ind = any(MeanEleValues~=0,2) & any(MeanEleValues==0,2);
    MeanEleValues(ind,:) = mean(MeanEleValues(ind,:)~=0);
    
end
