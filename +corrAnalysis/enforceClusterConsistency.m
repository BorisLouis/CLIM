function [ind2Add] = enforceClusterConsistency(ind2Add,currCluster,data,thresh)
    %function to enforce consistency within the cluster by checking that
    %to-be-added element are correlated with a random subset of the cluster
    nInd = length(ind2Add);
    nClust = length(currCluster);
    newList = ind2Add;
    if ~isempty(currCluster)
        
        %reshape input to compare correlation
%         id2Add  = find(ismember(indsCopy,newList));
% %         newList = repmat(id2Add,1,length(currCluster))';
% %         newList = newList(:);
%         currCluster = find(ismember(indsCopy,currCluster));
       % currCluster = repmat(currCluster,nInd,1);
        
        %get distance Map between new pixel and current cluster
        list2Comp = [newList;currCluster];
        
        distanceMap = corrAnalysis.getDistanceMapFromPxList(list2Comp,data);
        
        corrVal = distanceMap(nInd+1:end,1:nInd);
        
%         %combine input to extract value from distance map
%         subs = [newList,currCluster];    
%         inds = sub2ind(size(distanceMap),subs(:,1),subs(:,2));
%         corrVal = distanceMap(inds);
%         %reshape corrVal to test for each cluster
%         corrVal = reshape(corrVal,nClust,nInd);
        %check threshold
        corrTest = corrVal<thresh;
        %we want new data point to be correlated to 70% of the subset from the
        %existing group
        idx2delete = sum(corrTest,1)./size(corrTest,1)<0.7;

        %delete the data not matching the requested threshold
        ind2Add(idx2delete) = [];    
        %we dont want any duplicate of index (although it should
        %not happen)
        ind2Add = unique(ind2Add);

    end
end