function [corrMask,cleanCorrMask] = corrClustering(listCorrPx,sumPx,inds,data,thresh)
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
        %get Index of most correlated pixel
        [~,idx] = max(sumPx);
        %take the index of the first pixel to be treated
        currIndex = inds(idx);
        %add to a new list the pixel that are correlated with the
        %currently treated pixel
        currList = listCorrPx{idx};
        %check that the list does not contain already treated
        %pixels
        currList(ismember(currList,treatedIdx,'row'),:) =[];
        %if the list is empty then the pixel is "dead" and we
        %remove it from the list
        if isempty(currList)
            listCorrPx(idx) = [];
            inds(idx) = [];
            sumPx(idx) = [];
            idx = find(isnan(treatedIdx(:,1)),1);
            treatedIdx(idx) = currIndex;
            
        else
            %otherwise we treat all the pixel in the list in the
            %same way
            while ~isempty(currList)
                % add the new pixel to the group and keep track of it
                corrMask(currIndex) = group;
                clusters{group} = [clusters{group}; currIndex, find(indsCopy==currIndex)]; 
            
                idx = find(isnan(treatedIdx(:,1)),1);
                treatedIdx(idx) = currIndex;

                % remove it from the list
                listCorrPx(inds==currIndex) = [];
                sumPx(inds==currIndex)= [];
                inds(inds==currIndex) = [];
               
                %update the currentlist to add the new data
                currIndex = currList(1);
                currList(1,:) = [];
                
                %get list of element correlated to the current element
                list2Add  = listCorrPx{inds==currIndex}; 
                %remove already treated cases and cases already in the list
                list2Add(or(ismember(list2Add,treatedIdx,'row'),ismember(list2Add,currList,'row'))) = []; 
                
                if isempty(list2Add)
                    
                else
                    %find pixel that are already inside the cluster to check
                    %for correlation
                    currCluster = find(corrMask==group);

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
                    
                    %add the new element to the list
                    currList  = [currList;list2Add];
                    
                end

                if isempty(currList)
                    disp(['Found ' num2str(group) ' group(s).']);
                end

                count = count+1;
            end
            %mark the last case as treated:
            corrMask(currIndex) = group;
            idx = find(isnan(treatedIdx(:,1)),1);
            treatedIdx(idx) = currIndex;

            % remove it from the list
            listCorrPx(inds==currIndex) = [];
            sumPx(inds==currIndex) =[];
            inds(inds==currIndex) = [];

            %if exit the first while loop, the first group is complete, we need to
            %increment the group number as we are supposed to have treated all
            %cases.
            group = group+1;
            clusters{group} = [];
        end

        count=count+1;

        if count >= safeCount
            error('something went wrong')
        end
    end
    
    clusters(cellfun(@isempty,clusters)) = [];
    
    %fast cleanUp
    cleanCorrMask = rand(size(data,1),size(data,2));
%     try
%         distanceMap = corrAnalysis.getDistanceMapFromPxList(indsCopy,data);
%         
%     catch
%         
%         distanceMap = [];
%         warning('Too many pixels to generate distance map, using slower but more memory efficient route');
%     end
%     
%     if isempty(distanceMap)
%          [cleanCorrMask] = corrAnalysis.clusterCleanUpMemEff(corrMask,clusters,data,thresh);
%         
%     else
%         [cleanCorrMask] = corrAnalysis.clusterCleanUp(corrMask,clusters,distanceMap,thresh);
%        
%     end
            
end