function [clustEval,relNum] = evalClusters(mask,data)
    assert(ndims(data)==3,'Data is expected to have 3 dimensions');

    if iscell(mask)
        
        
    elseif ismatrix(mask)
        
        nClusters = max(mask(:));
        
        clustEval = table(zeros(nClusters,1),zeros(nClusters,1),zeros(nClusters,1),...
            zeros(nClusters,1), zeros(nClusters,1),'VariableNames',...
            {'meanCorr','medCorr','minCorr','maxCorr','nPixels'});
        
        for i = 1 :nClusters
           
            %get indices of mask for current cluster index
            currInds = find(mask==i);
            %if only one pixel we just give one's everywhere
            if length(currInds) ==1
                
                clustEval.meanCorr(i) = 1;
                clustEval.medCorr(i) = 1;
                clustEval.minCorr(i) = 1;
                clustEval.maxCorr(i) = 1;
                clustEval.nPixels(i) = 1;
                
            else
            %get Distance map for current cluster
            corrMat = 1-corrAnalysis.getDistanceMapFromPxList(currInds,data);
            
            %calculate some simple stats on the cluster
            clustEval.meanCorr(i) = mean(corrMat(corrMat<1));
            clustEval.medCorr(i) = median(corrMat(corrMat<1));
            clustEval.minCorr(i) = min(corrMat(corrMat<1));
            clustEval.maxCorr(i) = max(corrMat(corrMat<1));
            clustEval.nPixels(i) = length(currInds);
            
            
            end            
        end
        
        relNum.meanCorr = wmean(clustEval.meanCorr,clustEval.nPixels);
        relNum.medCorr  = wmean(clustEval.medCorr,clustEval.nPixels);
        relNum.minCorr  = wmean(clustEval.minCorr,clustEval.nPixels);
        relNum.maxCorr  = wmean(clustEval.maxCorr,clustEval.nPixels);
        relNum.nClusters = nClusters;
        
        
    else
        error('Unexpected mask format');
    end


end

function [wAvg] = wmean(data,weight)
    assert(isvector(data),'data is expected to be a vector');
    assert(isvector(weight),'data is expected to be a vector');
    assert(all(size(data) == size(weight)),'data and weight are expected to be the same size');
    
    wAvg = sum(data.*weight)/sum(weight);
    
end


