function [corrMask] = clusterCleanUp(corrMask,clusters,inds,distanceMap)

    nCluster = max(corrMask(:));
  
    for i = 1:nCluster
        currCluster = clusters{i};
        nPx = length(currCluster);
        
        for j = 1:nPx
            currentPx = currCluster(j,1);
            idx2Px    = currCluster(j,2);
            pxVec = repmat(idx2Px,nPx,1);
            
            idx   = sub2ind(size(distanceMap),pxVec,currCluster(:,2));
            correlation2Cluster = mean(distanceMap(idx));
            
            
            
        end
        
        
        
        
        
    end







end