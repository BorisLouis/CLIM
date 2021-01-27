function [corrMask] = corrClusteringV2(dim,inds,distanceMap)
    
    corrMask = zeros(dim(1),dim(2));
    clustCell = cell(1,100);
    delta = 1/44;
    
    pxList = inds;
    clustNum = 1;
    count = 0;
    while ~isempty(pxList)
        %Take a point
        currPoint = pxList(1);
        
        id = inds==currPoint;
        %get all the pixel correlated to that point and put them in a
        %cluster
        idx = distanceMap(id,:)==1;
        
        currentCluster = inds(idx);
        
        clusterList = find(idx);
        %check for "bad" pixels in the cluster
        for i = 1:length(clusterList)
            
            currentIdx = clusterList(i);
            
            positives = distanceMap(currentIdx,:)==1;
            positivesInCluster = and(positives,idx);
            positivesOutOfCluster = and(positives,~idx);
            
            %need to understand why 3 delta bad and 7 delta good, something
            %is off as a point can be both 3 delta bad and 7 delta good.
            deltaCriterion1 = or(~(sum(positivesInCluster) >= (1-3*delta)*sum(idx)),...
                    ~(sum(positivesOutOfCluster) <= 3*delta*sum(idx)));
            %TODO make criterion based on correlation level
            if deltaCriterion1
                
                %remove the point from the cluster
                currentCluster(currentCluster==inds(currentIdx)) = [];
           
            end
        end
        
        %check for good pixels outside the cluster
        [Y,idxA] = setdiff(pxList,currentCluster);
        
        idxB = zeros(size(idx));
        idxB(idxA)=1;
        for i = 1: length(Y)
            currentIdx = idxA(i);
            positives = distanceMap(currentIdx,:)==1;
            positivesInCluster = and(positives,idxB);
            positivesOutOfCluster = and(positives,~idxB);
            
            deltaCriterion2 = and((sum(positivesInCluster) >= (1-7*delta)*sum(idx)),...
                    (sum(positivesOutOfCluster) <= 7*delta*sum(idx))); 
            
            if deltaCriterion2
                
                %remove the point from the cluster
                currentCluster = [currentCluster; Y(i)];

            end
            
        end
        %add the cluster to the mask
        corrMask(currentCluster) = clustNum;
        
        clustCell{clustNum} = currentCluster;
        
        %update the cluster number
        clustNum = clustNum+1;
        
        %remove pixel of current cluster from pxList
        pxList(ismember(pxList,currentCluster)) = [];
        
        count = count +1;
        
        if count >length(inds)
            error('Something went wrong');
        end
        
    end
    

end