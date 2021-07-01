function [corrMask,allCorrMask] = bottomUpCorrClustering(corrRelation,data)
    %BottomUpClustering should provide the data and the threshold to
    %perform our standard correlation clustering
    
    corrStatus = Core.CorrStatus(corrRelation.listPx,corrRelation.listVal,...
        corrRelation.sumPx,corrRelation.indPx);
    
    minThresh  = corrRelation.tRange(1);
    safeCount = 2*size(data,1)*size(data,2);
    %initialize the variable
    count = 1;
    group = 1;
    corrMask = zeros(size(data,1),size(data,2));

    clusters{group} = [];
    
    %as long as the list of pixel that are correlated is not empty
    %we keep going
    while ~corrStatus.isFinished()
        %get Index of most correlated pixel
        [~,idx2MaxCorr] = min(corrStatus.sumPx);
        %take the index of the first pixel to be treated
        currIndex = corrStatus.inds(idx2MaxCorr);
        %add to a new list the pixels that are correlated with the
        %currently treated pixel
        refList = corrStatus.listPx{idx2MaxCorr};
        
        refVal  = corrStatus.listVal{idx2MaxCorr};
        %remove pixel that would have lower correlation than minimum threshold
        idx = refVal<minThresh;
        refList = refList(idx);
        refVal  = refVal(idx);
        
        %add these value to the cluster
        currentCluster = Core.Cluster(refList,refVal);
        
        %make the current List of pixel to be treated from the pixel
        %just added
        idx = ismember(corrStatus.inds,refList);
        currList = unique(cell2mat(corrStatus.listPx(idx)));
        
        %Delete added pixel from the list
        corrStatus.deletePxFromList(idx);        
        
        %check that the list does not contain already treated
        %pixels
        currList(ismember(currList,corrStatus.treatedPx,'row'),:) =[];
        
        if isempty(currList)
            listCorrPx(idx2MaxCorr) = [];
            inds(idx2MaxCorr) = [];
            sumPx(idx2MaxCorr) = [];
            idx = find(isnan(treatedIdx(:,1)),1);
            treatedIdx(idx) = currIndex;
            
        else
            %otherwise we treat all the pixel in the list in the
            %same way
            while ~isempty(currList)
                
                idx2Clust = a;
                
                
                
                
                
            end
        end
        
    
    
    
    end
            
end