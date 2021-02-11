function [corrMask,cleanCorrMask] = corrClustering(listCorrPx,inds,distanceMap,dim,thresh)
    %function actually does the clustering.
    safeCount = 2*dim(1)*dim(2);
    %initialize the variable
    count = 1;
    group = 1;
    corrMask = zeros(dim(1),dim(2));
    treatedIdx = NaN(length(inds),1);
    indsCopy = inds;
    indsPxIdx = 1:length(inds);
    clusters{group} = [];
    %as long as the list of pixel that are correlated is not empty
    %we keep going
    while ~isempty(listCorrPx) 
        %take the index of the first pixel to be treated
        currIndex = inds(1);
        %add to a new list the pixel that are correlated with the
        %currently treated pixel
        currList = listCorrPx{1};
        %check that the list does not contain already treated
        %pixels
        currList(ismember(currList,treatedIdx,'row'),:) =[];
        %if the list is empty then the pixel is "dead" and we
        %remove it from the list
        if isempty(currList)
            listCorrPx(1) = [];
            idx = find(isnan(treatedIdx(:,1)),1);
            treatedIdx(idx) = currIndex;
            inds(1) = [];
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
                inds(inds==currIndex) = [];

                %update the currentlist to add the new data
                currIndex = currList(1);
                currList(1,:) = [];
                %Quick bug fix is to check if currIndex is in inds
                %and if not we just delete it?

                %get list of element correlated to the current element
                list2Add  = listCorrPx{inds==currIndex}; 
                %remove already treated cases and cases already in the list
                list2Add(or(ismember(list2Add,treatedIdx,'row'),ismember(list2Add,currList,'row'))) = []; 

                %find pixel that are already inside the cluster to check
                %for correlation
                currCluster = find(corrMask==group);

                %take a random subset of 10 pixel from the
                %current cluster
                try
                    r = randperm(length(currCluster,10));
                    currCluster = currCluster(r);
                catch

                end
                %check that the pixel that should be added to the
                %list are indeed correlated to the current cluster
                [list2Add] = corrAnalysis.enforceClusterConsistency(list2Add,...
                    currCluster,distanceMap,indsCopy,thresh);
                %add the new element to the list
                currList  = [currList;list2Add];

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
    %test loop clean up
    
    [cleanCorrMask,clusters] = corrAnalysis.clusterCleanUp(corrMask,clusters,distanceMap);
    
%     for i = 1:10
%         [cleanCorrMask,clusters] = corrAnalysis.clusterCleanUp(cleanCorrMask,clusters,distanceMap);
%         
%     end
    
    
    %[cleanCorrMask] = corrAnalysis.clusterCleanUp2(corrMask,clusters,distanceMap);
    
end