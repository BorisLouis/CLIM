function [clustEval,relNum,corr2ClusterMap] = evalClusters(mask,data)
    assert(ndims(data)==3,'Data is expected to have 3 dimensions');

    if iscell(mask)
        
        
    elseif ismatrix(mask)
        
        clusterID = unique(mask(mask>0));
        nClusters = length(clusterID);
        
        clustEval = table(zeros(nClusters,1),zeros(nClusters,1),zeros(nClusters,1),...
            zeros(nClusters,1), zeros(nClusters,1),zeros(nClusters,1),...
            zeros(nClusters,1),'VariableNames',...
            {'meanCorr','medCorr','std','varCoeff','minCorr','maxCorr','nPixels'});
        traces = zeros(nClusters,size(data,3));
        corr2ClusterMap = zeros(size(mask));
        for i = 1 :nClusters
           
            %get indices of mask for current cluster index
            currInds = find(mask==clusterID(i));
            %if only one pixel we just give one's everywhere
            if length(currInds) ==1
                
                clustEval.meanCorr(i) = 1;
                clustEval.medCorr(i) = 1;
                clustEval.minCorr(i) = 1;
                clustEval.maxCorr(i) = 1;
                clustEval.nPixels(i) = 1;
                
            else
                % We calculate the average trace
                [row,col] = find(mask==clusterID(i));
                r = repmat(row,1,size(data,3));
                r = reshape(r',size(data,3)*length(row),1);
                c = repmat(col,1,size(data,3));
                c = reshape(c',size(data,3)*length(col),1);
                
                f = repmat((1:size(data,3))',length(col),1);
                idx = sub2ind(size(data),r,c,f);
                
                tmpTrace = data(idx);
                tmpTrace = reshape(tmpTrace,size(data,3),length(row));
                
                traces(i,:) = mean(tmpTrace,2);
                
                %get Distance map for current cluster
                corrMat = 1-corrAnalysis.getDistanceMapFromPxList(currInds,data);
                
                corrMat(corrMat==1) = nan;
                
                corr2ClusterMap(currInds) = nanmean(corrMat,2);
                
                %calculate some simple stats on the cluster
                clustEval.meanCorr(i) = nanmean(corrMat(corrMat<1));
                clustEval.medCorr(i)  = nanmedian(corrMat(corrMat<1));
                clustEval.std(i)      = nanstd(corrMat(corrMat<1));
                clustEval.varCoeff(i) = clustEval.std(i)/clustEval.meanCorr(i);
                clustEval.minCorr(i) = nanmin(corrMat(corrMat<1));
                clustEval.maxCorr(i) = nanmax(corrMat(corrMat<1));
                clustEval.Intensity(i) = mean(traces(i,:));
                clustEval.nPixels(i) = length(currInds);
            
            end            
        end
        %Process the traces
        corrCoeff = corrcoef(traces');
        
        U = triu(corrCoeff);
        U = nonzeros(U);
        data = U(U<1);
        
        % store interclusterCorr
        corrCoeff(corrCoeff==1) = nan;
        clustEval.meanInterClustCorr = nanmean(corrCoeff,2);
        clustEval.maxInterClustCorr  = nanmax(corrCoeff,[],2);
        
        
        relNum.meanCorr = wmean(clustEval.meanCorr,clustEval.nPixels);
        relNum.medCorr  = wmean(clustEval.medCorr,clustEval.nPixels);
        relNum.std      = wmean(clustEval.std,clustEval.nPixels);
        relNum.varCoeff = wmean(clustEval.varCoeff,clustEval.nPixels);
        relNum.minCorr  = wmean(clustEval.minCorr,clustEval.nPixels);
        relNum.maxCorr  = wmean(clustEval.maxCorr,clustEval.nPixels);
        relNum.meanInt  = wmean(clustEval.Intensity,clustEval.nPixels);
        relNum.meanSize = mean(clustEval.nPixels);
        relNum.meanInterClusterCorr = mean(data);
        relNum.stdInterClusterCorr  = std(data);
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

