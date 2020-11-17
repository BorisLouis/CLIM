function [cleanCorrMask] = cleanCorrMask(data, corrMask, signThresh)
    nClusters = max(corrMask(:));
    nFrames = size(data,3);
    allTraces = cell(2,nClusters);
    meanTraces = zeros(nClusters,nFrames);
    
    for i = 1: nClusters
       BW = corrMask==i;
       BW = bwareaopen(BW,4,8);
       mBW = repmat(BW,1,1,nFrames);
       trace = data(mBW);
       traces = reshape(trace',sum(BW(:)),nFrames);
       allTraces{1,i} = traces; 
       allTraces{2,i} = mean(traces,1);
       
       meanTraces(i,:) = allTraces{2,i};
       
    end
    
    for i = 1:nClusters
       
        currTraces = allTraces{1,i};
        %stack the two together to calculate correlation
        currTraces = [currTraces; meanTraces];
        corr = corrcoef(currTraces');
        
        corr2Compare = corr(end-(nClusters-1):end,:);
        
          
        idx = corr2Compare>signThresh*corr2Compare(i,:);
        test = sum(idx,1);
        id   = find(test>1);
        %delete excess idx
        id(id>size(allTraces{1,i},1)) = [];
        
        trace2Del = allTraces{1,i}(id,1);
        idx2Del = ismember(data(:,:,1),trace2Del);
        
        corrMask(idx2Del) = 0;
        
        %need to figure out the link between the order of the trace and
        %corrMask so we can delete the correct indices
        
        
    end
    
    cleanCorrMask = zeros(size(corrMask));
    for i = 1:nClusters
       
        BWCopy = corrMask==i;
        
        BWCopy = bwareaopen(BWCopy,4,8);
        
        cleanCorrMask(BWCopy) = max(cleanCorrMask(:))+1;
        
    end



end