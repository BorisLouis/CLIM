function [corrMask] = corrClusteringBestOfThree(meanPx,inds,corrMat,dim)
    %The idea of the modifications here is to:
    %1) Use listVal to and threshold to add to cluster(this will allow to
    %test various threshold without having to re-run the first steps)
    %2) Add multiple indices to the cluster at once instead of 1 by 1 to
    %reduce the processing time.
    kThresh = 3;
    %function actually does the clustering.
    safeCount = 2*size(corrMat,1);
    %initialize the variable
    count = 1;
    group = 1;
    corrMask = zeros(dim(1),dim(2));
    treatedIdx = NaN(length(inds),1);
    indsCopy = inds;
    
    clusters{group} = [];
    %as long as the list of pixel that are correlated is not empty
    %we keep going
    while ~all(meanPx == 0) 
        
        %get Index of most correlated pixel m
        [~,idx] = max(meanPx);
        %take the index of the first pixel to be treated
        currIndex = inds(idx);
        %add to a new list the pixel that are correlated with the
        %currently treated pixel
        currCorrList = corrMat(idx,:);
        
        %get highest correlated pixels:
        maxK = maxk(currCorrList,kThresh+1);
        
       
        currList = inds(currCorrList>=min(maxK));
        
        %check that the list does not contain already treated
        %pixels
        %currVal(ismember(currList,treatedIdx,'row'),:) =[];
        currList(ismember(currList,treatedIdx,'row'),:) =[];

        %add the pixel to the cluster
        %corrMask(tmpList) = group;
        
        clusters{group} = [clusters{group}; currList];
        %keep track of the added pixels
        idx = find(isnan(treatedIdx(:,1)),1);
        treatedIdx(idx:idx+length(currList)-1) = currList;
        
        %if the list is empty then the pixel is "dead" and we
        %remove it from the list
        if isempty(currList) 
            meanPx(ismember(inds,currList)) = 0;
            idx = find(isnan(treatedIdx(:,1)),1);
            treatedIdx(idx:idx+length(currList)-1) = currList;
            
        else
            %otherwise we treat all the pixel in the list in the
            %same way
            while ~isempty(currList)
                
                %1) Treat the list (get correlated pixels and add to treated pixels)
                currCorrList = corrMat(sort(currList),:);
                currMaxK = maxk(currCorrList,kThresh+1,2);
                idx = sum(currCorrList>=min(currMaxK,[],2),1);
                newCurrList = inds(idx>0);

                % remove it from the list
                meanPx(ismember(inds,currList)) = 0;

                %2) Renew the list
                currList = newCurrList;     
                
                %3) remove already treated Idx:
                currList(ismember(currList,treatedIdx,'row')) = [];
                
                if isempty(currList)
                    disp(['Found ' num2str(group) ' group(s).']);
                else
                    %mark as treated index
                    idx = find(isnan(treatedIdx(:,1)),1);
                    treatedIdx(idx:idx+length(currList)-1) = currList;
                end

                count = count+1;
                
                %add the pixel to the cluster
                %corrMask(currList) = group;
                clusters{group} = [clusters{group}; currList];
                
               
            end
            %if exit the first while loop, the first group is complete, we need to
            %increment the group number as we are supposed to have treated all
            %cases.
            
            if length(clusters{group})<5
                clusters{group} = [];
            else
                corrMask(clusters{group}(:,1)) = group;
                group = group+1;
                clusters{group} = [];
            end

            count=count+1;

            if count >= safeCount
                error('something went wrong')
            end
        end                
    end