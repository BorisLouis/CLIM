function [corrMask] = corrClusteringV4(data,corrRelation,distanceMap,thresh)
    inds = corrRelation.indPx;
    dim = size(data);
    clusterList = Core.CorrClusterList(data,inds);
    clusterList.updateState
    distMap = distanceMap;
    clear data;
    
    nCounts = 1;
    stopCriteria = false;   
    
    while not(stopCriteria) 

            %get Index of most correlated pixel m
            [val,idx] = min(distanceMap(distanceMap>0));
            id = find(distanceMap==val);
            [i,j] = ind2sub(size(distanceMap),id(1));
            
            %merge the two clusters i,j
            [didMerge] = clusterList.merge(i,j,distMap);
            
            if not(didMerge)
                nCounts = nCounts + 1;
                distanceMap(i,j) = inf;
                distanceMap(j,i) = inf;
                
            else
                %modify distance map to take into account the new cluster
                distMapA = distanceMap(i,:);
                distMapB = distanceMap(j,:);
                
                %!!! not weighted at the moment
                meanMap = (distMapA+distMapB)/2;
                
                distanceMap(i,:) = [];
                distanceMap(j,:) = [];
                distanceMap(:,i) = [];
                distanceMap(:,j) = [];
                
                meanMap(i) = [];
                meanMap(j) = [];
                
                distanceMap(end+1,:) = meanMap;
                meanMap = [meanMap 0];
                distanceMap(:,end+1) = meanMap;
                
                nCounts = 1; 
            end
            
            stopCriteria = or(nCounts> length(clusterList.clusters),length(distanceMap)==1);
    end

    corrMask = clusterList.getCorrMask(dim);


end
