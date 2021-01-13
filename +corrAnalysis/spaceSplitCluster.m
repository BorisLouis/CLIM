function [hierarchical] = spaceSplitCluster(corrMask)
    %Split cluster that have the same number but are spatially separated
    %into two different subclusters.
    
    hierarchical = cell(size(corrMask));
    for i = 1:max(corrMask(:))
         
        %get the requested cluster
        subMask = corrMask == i;
        currMask = bwlabel(subMask);

        idx = max(currMask(:));
        if idx == 1
            %there is no cluster separated spatially, no need to split
            hierarchical(subMask>0) = [];
                        
        else
            %we loop through the different separated cluster to give them
            %an additional index
            
            for j = 1:idx
                subMask = currMask==j;
                
                %Keep track of clusters hierarchy
                idxMat = [ones(length(subMask(subMask>0)),1)*j];

                hierarchical(subMask>0) = num2cell(idxMat);

            end  
        end
            
    end
end