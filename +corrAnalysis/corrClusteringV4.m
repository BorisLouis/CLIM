function [corrMask] = corrClusteringV4(data,corrRelation,distanceMap,thresh)
    inds = corrRelation.indPx;
        
    clusterList = Core.CorrClusterList(data,inds);
    clusterList.updateState
    
    clear data;
    
    nCounts = 1;
    stopCriteria = false;   
    
    while not(stopCriteria) 

            %get Index of most correlated pixel m
            [val,idx] = min(distanceMap(distanceMap>0));
            
            [i,j] = ind2sub(size(distanceMap),idx(1));
            
            
            clusterList.merge(i,j,distanceMap);
            
            
            %take the index of the first pixel to be treated
            currIndex = inds(i);
            %add to a new list the pixel that are correlated with the
            %currently treated pixel
            currList = listCorrPx{idx};
            currVal  = listVal{idx};
            
            
            
            
            nCount
    end



end
