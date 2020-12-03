function [hierarchical] = spaceSplitCluster(corrMask)
    hierarchical = cell(size(corrMask));
    for i = 1:max(corrMask(:))
         
        %get the requested cluster
        subMask = corrMask == i;
        currMask = bwlabel(subMask);

        idx = max(currMask(:));

        for j = 1:idx
            subMask = currMask==j;
            %store the cleaned mask
            corrMask(corrMask==i) =0;
            corrMask(subMask>0) = i;
            %Keep track of clusters hierarchy
            nElem = ones(length(subMask(subMask>0)),1);
            idxMat = [nElem*i,nElem*j];
           
            hierarchical(subMask>0) = mat2cell(idxMat,nElem);
            
        end   
    end
end