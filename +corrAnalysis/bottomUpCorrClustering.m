function [corrMap,clusters] = bottomUpCorrClustering(corrRelation,data)
    %BottomUpClustering should provide the data and the threshold to
    %perform our standard correlation clustering
    
    corrStatus = Core.CorrStatus(corrRelation.listPx,corrRelation.listVal,...
        corrRelation.sumPx,corrRelation.indPx);
    
    minThresh  = corrRelation.tRange(1);
    safeCount = 2*size(data,1)*size(data,2);
    %initialize the variable
    count = 1;
    group = 1;
    corrMap = zeros(size(data,1),size(data,2));

    clusters{group} = [];
    
    %as long as the list of pixel that are correlated is not empty
    %we keep going
    while ~corrStatus.isFinished()
        if count>safeCount
            error('something went wrong')
        end
        %get Index of most correlated pixel
        [~,idx2MaxCorr] = min(corrStatus.sumPx);
        %get the corresponding index in the image space
        currIndex = corrStatus.inds(idx2MaxCorr);
        currVal   = corrStatus.listVal{idx2MaxCorr};
        
        %create a new cluster base on that pixel
        currentCluster = Core.Cluster(currIndex,0,group);
        
        %add to a new list the pixels that are correlated with the
        %currently treated pixel(s)
        currList = corrStatus.listPx{idx2MaxCorr};
        
        currListVal  = corrStatus.listVal{idx2MaxCorr};
                   
        %remove pixel that would have lower correlation than minimum threshold
        idx = currListVal<minThresh;
        currList = currList(idx);
        currListVal  = currListVal(idx);
        
        %remove the current pixel from the list of available pixel
        idx = ismember(corrStatus.inds,currIndex);
        corrStatus.deletePxFromList(idx);
        
        %check that the list does not contain already treated
        %pixels
        currListVal(ismember(currList,corrStatus.treatedPx,'row'),:) = [];
        currList(ismember(currList,corrStatus.treatedPx,'row'),:) = [];
         
        %otherwise we treat all the pixel in the list in the
        %same way
        while ~isempty(currList)
            
            %add these value to the cluster
            listOfAddedPx = currentCluster.addPxToList(currList,currListVal);

            %make the current List of pixel to be treated from the pixel
            %just added
            idx = ismember(corrStatus.inds,listOfAddedPx);
            currList = unique(cell2mat(corrStatus.listPx(idx)));
            %check that the list does not contain already treated
            %pixels
            currList(ismember(currList,corrStatus.treatedPx,'row'),:) =[];
            
            [currListVal] = getCorrVal(currList,currentCluster.pxList(1),data);
                        
            %Delete added pixel from the list
            corrStatus.deletePxFromList(idx);
            
            %remove pixel that would have lower correlation than minimum threshold
            idx = currListVal<minThresh;
            currList = currList(idx);
            currListVal  = currListVal(idx);
                        
            
               
            

            count = count+1;
            
        end
        %if exit the first while loop, the current group is complete, we need to
        %increment the group number as we are supposed to have treated all
        %cases.
        %store the current cluster
        if length(currentCluster.pxList)>10
            currentCluster.print;
            clusters{group} = currentCluster;
            corrMap(currentCluster.pxList) = group;
            group = group+1;
            clusters{group} = [];
        end
    
    end
end

function [corrVal] = getCorrVal(pxList,refPx,data)
    
    %append the list of pixel to the reference pixel
    px2Corr = [refPx;pxList];
    
    %extract traces from the list of px
    [row,col] = ind2sub(size(data),px2Corr);
    r = repmat(row,1,size(data,3));
    r = reshape(r',size(data,3)*length(row),1);
    c = repmat(col,1,size(data,3));
    c = reshape(c',size(data,3)*length(col),1);

    f = repmat((1:size(data,3))',length(col),1);
    idx = sub2ind(size(data),r,c,f);

    %get the data
    tmpTrace = data(idx);
    tmpTrace = double(reshape(tmpTrace,size(data,3),length(row)));
    
    corrVal = 1-corr(tmpTrace(:,1),tmpTrace(:,2:end));
    
    corrVal = corrVal(:);
    
    
end