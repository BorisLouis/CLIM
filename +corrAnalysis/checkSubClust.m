function [checkRes] = checkSubClust(MLCorrMask,data,stopCriteria)
    
    nSubClust = max(MLCorrMask(:));
    %calculate main cluster distanceMap
    mainClust = MLCorrMask>0;
    mainInds = find(mainClust);
    mainDistMap = corrAnalysis.getDistanceMapFromPxList(mainInds,data);
    
    variance = zeros(nSubClust,1);
    stopMetric = variance;
    
    variance(1) = var(mainDistMap);
    stopMetric(1) = sqrt(sum(mainClust(:)))/variance(1);
    
    for i = 1:nSubClust
        tmpClust = MLCorrMask == i;
        tmpInds = find(tmpClust);
        tmpDistMap = corrAnalysis.getDistanceMapFromPxList(tmpInds,data);

        variance(i+1) = var(tmpDistMap);
        stopMatric(i+1) = sqrt(sum(tmpClust(:)))/variance(i+1);
        
        
    end
    
    
    %test
    


end