function [checkRes] = checkSubClust(MLCorrMask,data)
    
    nSubClust = max(MLCorrMask(:));
    %calculate main cluster distanceMap
    mainClust = MLCorrMask>0;
    mainInds = find(mainClust);
    mainDistMap = corrAnalysis.getDistanceMapFromPxList(mainInds,data);
    
    
    %Criteria 1 ==> CV needs to decrease (smaller STD and/or larger mean)
    stdev = zeros(nSubClust+1,1);
    CV = stdev;
    stdev(1) = std(mainDistMap(:));
    CV(1) = stdev(1)/(1-mean(mainDistMap(:)));
    
    for i = 1:nSubClust
        %Loop through the sub cluster
        tmpClust = double(MLCorrMask == i);        
        tmpInds = find(tmpClust);
        %get the distance map
        tmpDistMap = corrAnalysis.getDistanceMapFromPxList(tmpInds,data);
        
        %Calculate the coefficient of variation
        stdev(i+1) = std(tmpDistMap(:));
        CV(i+1) = stdev(i+1)/(1-mean(tmpDistMap(:)));
        
        
    end
    %Here we check if the coefficient of variation improves for all sub
    %cluster
    
    if all(CV(2:3) < CV(1))
        
        checkRes = true;
    
    else
        
        checkRes = false;
        
    end
    
end