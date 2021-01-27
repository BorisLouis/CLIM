function [corrMask] = corrClusteringV3(dim,inds,distanceMap,corrThreshold)
    
    corrMask = zeros(dim(1),dim(2));
    clustCell = cell(1,100);
    pxList = 1:length(inds);
    clustNum = 1;
    count = 0;
    BWDistanceMap = distanceMap<corrThreshold;
    while ~isempty(pxList)
        %Take a point
        
        currMap = distanceMap(pxList,:);
        currMap = currMap(:,pxList);
        avgCorr = mean(currMap,1);
        [~,minAvgIdx] = min(avgCorr);
        
        %get all the pixel correlated to that point and put them in a
        %cluster
        currBWMap = BWDistanceMap(pxList,:);
        currBWMap = currBWMap(:,pxList);
        idx = currBWMap(minAvgIdx,:)==1;
        
        currentCluster = inds(pxList(idx));
       
        %add the cluster to the mask
        corrMask(currentCluster) = clustNum;
        
        clustCell{clustNum} = [currentCluster pxList(idx)'];
        
        %update the cluster number
        clustNum = clustNum+1;
        
        %remove pixel of current cluster from pxList
        pxList(ismember(inds(pxList),currentCluster)) = [];
        
        count = count +1;
        
        if count >length(inds)
            error('Something went wrong');
        end
        
    end
    
    [corrMask] = corrAnalysis.clusterCleanUp(corrMask,clustCell,inds,distanceMap);
    
    
end