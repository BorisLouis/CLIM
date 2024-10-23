function [corrMask,frames] = corrClustering(listCorrPx,listVal,meanPx,inds,data,thresh,doPlot)
    %The idea of the modifications here is to:
    %1) Use listVal to and threshold to add to cluster(this will allow to
    %test various threshold without having to re-run the first steps)
    %2) Add multiple indices to the cluster at once instead of 1 by 1 to
    %reduce the processing time.

    %function actually does the clustering.
    safeCount = 2*size(data,1)*size(data,2);
    %initialize the variable
    count = 1;
    group = 1;
    corrMask = zeros(size(data,1),size(data,2));
    treatedIdx = NaN(length(inds),1);
    indsCopy = inds;
    clusters{group} = [];
    %as long as the list of pixel that are correlated is not empty
    %we keep going    
    while ~isempty(listCorrPx) 
        
        %get Index of most correlated pixel m
        [~,idx] = max(meanPx);
        %take the index of the first pixel to be treated
        currIndex = inds(idx);
        %add to a new list the pixel that are correlated with the
        %currently treated pixel
        currList = listCorrPx{idx};
        currVal  = listVal{idx};
        
        %remove value that have been treated:
        currVal(or(ismember(currList,treatedIdx,'row'),~ismember(currList,indsCopy,'row')),:) =[];
     
        currList(or(ismember(currList,treatedIdx,'row'),~ismember(currList,indsCopy,'row')),:) =[];
        

        
        currList(currVal<thresh) = [];
        %check that the list does not contain already treated
        %pixels
       
        
        tmpList = [currIndex; currList];
        
        %add the pixel to the cluster
        %corrMask(tmpList) = group;
        
        clusters{group} = [clusters{group}; tmpList];
        %keep track of the added pixels
        idx = find(isnan(treatedIdx(:,1)),1);
        treatedIdx(idx:idx+length(tmpList)-1) = tmpList;
        
        currList = tmpList;
        %if the list is empty then the pixel is "dead" and we
        %remove it from the list
        if isempty(currList)
            listCorrPx(ismember(inds,currList)) = [];
            listVal(ismember(inds,currList)) = [];
            meanPx(ismember(inds,currList)) = [];
            inds(ismember(inds,currList)) = [];
            idx = find(isnan(treatedIdx(:,1)),1);
            treatedIdx(idx:idx+length(currList)-1) = currList;
            
        else
            %otherwise we treat all the pixel in the list in the
            %same way
            while ~isempty(currList)
                
                %1) Treat the list (get correlated pixels and add to treated pixels)
                newCurrList = cell2mat(listCorrPx(ismember(inds,currList)));
                valList     = cell2mat(listVal(ismember(inds,currList)));
                
                % remove it from the list
                listCorrPx(ismember(inds,currList)) = [];
                listVal(ismember(inds,currList)) = [];
                meanPx(ismember(inds,currList)) = [];
                inds(ismember(inds,currList)) = [];

                %2) Renew the list
                currList = newCurrList;
                
                
                %3) apply restrictions to the list
                % correlation needs to be high enough:
                currList(valList<thresh) = [];
                %All pixel should only be there once:
                list2Add = unique(currList);
                
                %remove already treated pixels from the list:
                if ~isempty(list2Add)
                    list2Add(ismember(list2Add,treatedIdx,'row')) = []; 
                    
                end
                %% STOPPPED HERE
                if isempty(list2Add)
                                     
                else
                    list2Add(~ismember(list2Add,indsCopy,'row'),:) = [];
                    if ~isempty(list2Add)
                        %find pixel that are already inside the cluster to check
                        %for correlation
                        currCluster = clusters{group}(:,1);
                        %currCluster = find(corrMask==group);

                        %take a random subset of 10 pixel from the
                        %current cluster
                        try
                            r = randperm(length(currCluster),10);
                            currCluster = currCluster(r);
                        catch

                        end
                        %check that the pixel that should be added to the
                        %list are indeed correlated to the current cluster
                        [list2Add] = corrAnalysis.enforceClusterConsistency(list2Add,...
                            currCluster,data,thresh);
                    end
                end
                
                %add the new element to the list
                currList  = list2Add;
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

            if and(count >= safeCount,~isempty(currList))
                error('something went wrong')
            end
        end                
    end